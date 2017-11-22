function [gain_dB,fc,bw_oct,Ho,f] = parametricEQest( initGain, initFc, initBw, fs, H_mic, startF, endF, f_interp)
%parametricEQest - Estimate the parameters of the filter
%
    e = sum(abs(H_mic));
    gain_dB = initGain;
    fc = initFc;
    bw_oct = initBw;
    
    % parametric filter with initial values
    [Ho,~]=parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
    % initial error
    if(sum(abs(H_mic+Ho)) < e)
        e = sum(abs(H_mic+Ho));
    end
    
    for j = 1:10
        % fc calibration: from start freq to end freq (uniform (continuous) distribution)
        for i = 1:20
            param = rand * (endF-startF) + startF;
            [Ho,~]=parametricEQ(gain_dB,param,bw_oct,fs,f_interp);
            if (sum(abs(H_mic+Ho)) < e)
                e = sum(abs(H_mic+Ho));
                fc = param;
            end
        end
        % gain calibration: form 0.15*MaxGain to 1.15*MaxGain (MaxGain = initGain)
        for i = 1:20
            param = rand + 0.15;
            [Ho,~]=parametricEQ(initGain*param,fc,bw_oct,fs,f_interp);
            if (sum(abs(H_mic+Ho)) < e)
                e = sum(abs(H_mic+Ho));
                gain_dB = initGain*param;
            end
        end
        % bw calibration:  form 0.2*initBw to 1.2*initBw (initBw = bandwidth of the error area)
        for i = 1:20
            param = rand + 0.2;
            [Ho,~]=parametricEQ(gain_dB,fc,initBw*param,fs,f_interp);
            if (sum(abs(H_mic+Ho)) < e)
                e = sum(abs(H_mic+Ho));
                bw_oct = initBw*param;
            end
        end
    end
    
    % final filter
    [Ho,f]=parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
   
end

