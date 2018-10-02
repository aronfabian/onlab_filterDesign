function [Hb1,Hb2,fcb1,fcb2] = butterworthFilters(H_mic,f_interp,tol_interp,H_trgt,configFile)
%butterworthFilters - design a HP and LP Butterworth Filter 
%
    load(configFile)
    % HPF:
    fs = 44100;
    fcb1 = START_FREQ; % highpass filter cutoff freq
    [b,a] = butter(2,fcb1/(fs/2), 'high');
    [Hb1,~] = freqz(b,a,f_interp,fs);
    Hb1 = 20*log10(abs(Hb1));
%     e = sum(abs(H_mic+Hb1));
    [eA,eT] = errorCalc(H_mic,H_trgt,tol_interp,configFile);
    [eA_new,eT_new] = errorCalc(H_mic+Hb1,H_trgt,tol_interp,configFile);
    if(eT_new < eT)
        if(eA_new < eA)
            eA = eA_new;
            eT = eT_new;
        end
    end
    
    freq = logspace(log10(START_FREQ),log10(END_FREQ),100);
    for i = 2:length(freq)
        f1 = freq(i); % highpass filter cutoff freq
        [b,a] = butter(2,f1/(fs/2), 'high');
        [H1,~] = freqz(b,a,f_interp,fs);
        H1 = 20*log10(abs(H1));
%         
%         if(sum(abs(H_mic+H1)) < e)
%             e = sum(abs(H_mic+H1));
%             Hb1 = H1;
%             fcb1 = f1;
%         end

        [eA_new,eT_new] = errorCalc(H_mic+H1,H_trgt,tol_interp,configFile);
        if(eT_new < eT)
            if(eA_new < eA)
                eA = eA_new;
                eT = eT_new;
                Hb1 = H1;
                fcb1 = f1;
            end
        end
        
    end
    
    % if error bigger with filter then filter = 0
%     if (e > sum(abs(H_mic)))
%         Hb1 = zeros(1,length(f_interp));
%         fcb1 = Inf;
%     end
    [eA_new,eT_new] = errorCalc(H_mic,H_trgt,tol_interp,configFile);
    if(eT_new < eT)
        if(eA_new < eA)
            eA = eA_new;
            eT = eT_new;
            Hb1 = zeros(1,length(f_interp));
            fcb1 = Inf;
        end
    end
    
    
    % LPF:
    fcb2 = END_FREQ;
    [b,a] = butter(2,fcb2/(fs/2), 'low');
    [Hb2,~] = freqz(b,a,f_interp,fs);
    Hb2 = 20*log10(abs(Hb2));
    
%     e = sum(abs(H_mic+Hb1+Hb2));
    [eA_new,eT_new] = errorCalc(H_mic+Hb1+Hb2,H_trgt,tol_interp,configFile);
    if(eT_new < eT)
        if(eA_new < eA)
            eA = eA_new;
            eT = eT_new;
        end
    end
    
    for i = 1:(length(freq)-1)
        f2 = freq(end-i); % lowpass filter cutoff freq
        [b,a] = butter(2,f2/(fs/2), 'low');
        [H2,~] = freqz(b,a,f_interp,fs);
        H2 = 20*log10(abs(H2));
%         if(sum(abs(H_mic+Hb1+H2)) < e)
%             e = sum(abs(H_mic+Hb1+H2));
%             Hb2 = H2;
%             fcb2 = f2;
%         end 
        [eA_new,eT_new] = errorCalc(H_mic+Hb1+H2,H_trgt,tol_interp,configFile);
        if(eT_new < eT)
            if(eA_new < eA)
                eA = eA_new;
                eT = eT_new;
                Hb2 = H2;
                fcb2 = f2;
            end
        end
    end
    
    % if error bigger with filter then filter = 0
%     if (e > sum(abs(H_mic+Hb1)))
%         Hb2 = zeros(1,length(f_interp));
%         fcb2 = Inf;
%     end
    [eA_new,eT_new] = errorCalc(H_mic,H_trgt,tol_interp,configFile);
    if(eT_new < eT)
        if(eA_new < eA)
            eA = eA_new;
            eT = eT_new;
            Hb2 = zeros(1,length(f_interp));
            fcb2 = Inf;
        end
    end
    
end