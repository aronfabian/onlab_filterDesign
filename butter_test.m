n = 2;
Wn = 0.05;
ftype = 'high';

% Transfer Function design
[b,a] = butter(n,Wn,ftype);      % This is an unstable filter


% Display and compare results
%hfvt = fvtool(b,a,'FrequencyScale','log');
%legend(hfvt,'TF Design','ZPK Design')
sys = tf(a,b)
hpfilter = abs(freqresp(sys,f_interp_plot));
hp=squeeze(hpfilter);
hp=20*log10(hp);
semilogx(f_interp_plot,hp)