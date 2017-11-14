function [Hb1,Hb2] = butterworthFilters(H_mic, f_interp,START_FREQ,END_FREQ)
%butterworthFilters - design a HP and LP Butterworth Filter 
%
    % HPF:
    fs = 44100;
    fcb1 = START_FREQ; % highpass filter cutoff freq
    [b,a] = butter(1,fcb1/(fs/2), 'high');
    [Hb1,f] = freqz(b,a,f_interp,fs);
    Hb1 = 20*log10(abs(Hb1));
    e = sum(abs(H_mic+Hb1));
    
    freq = logspace(log10(START_FREQ),log10(END_FREQ),100);
    for i = 2:length(freq)
        f1 = freq(i); % highpass filter cutoff freq
        [b,a] = butter(1,f1/(fs/2), 'high');
        [H1,f] = freqz(b,a,f_interp,fs);
        H1 = 20*log10(abs(H1));
        
        if(sum(abs(H_mic+H1)) < e)
            e = sum(abs(H_mic+H1));
            fcb1 = f1;
            Hb1 = H1;
        end
        
    end
    
    % LPF:
    fcb2 = END_FREQ;
    [b,a] = butter(1,fcb2/(fs/2), 'low');
    [Hb2,f] = freqz(b,a,f_interp,fs);
    Hb2 = 20*log10(abs(Hb2));
    
    e = sum(abs(H_mic+Hb1+Hb2));
    
    for i = 1:(length(freq)-1)
        f2 = freq(end-i); % lowpass filter cutoff freq
        [b,a] = butter(1,f2/(fs/2), 'low');
        [H2,f] = freqz(b,a,f_interp,fs);
        H2 = 20*log10(abs(H2));
        
        if(sum(abs(H_mic+Hb1+H2)) < e)
            e = sum(abs(H_mic+Hb1+H2));
            fcb2 = f2;
            Hb2 = H2;
        end 
    end
end