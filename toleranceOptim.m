function [gain_dB,fc,bw_oct,Ho,f] = toleranceOptim( initGain, initFc, initBw, fs, H_mic, f_interp,COMP_FLINES, tol_interp, configFile)
%toleranceOptim - modifify the parameters of the parametric filter
%therefore the tansfer function is in the middle of the tolerance band
%
    load(configFile)
    
    gain_dB = initGain;
    fc = initFc;
    bw_oct = initBw;

    % effective range of the filter
%     limit = 0.1;
    [H_filt,f] = parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
%     rangeHigh = find(H_filt(1:end-1)>limit & H_filt(2:end)<limit);
%     rangeLow = find(H_filt(1:end-1)<limit & H_filt(2:end)>limit);
%     

    rangeHigh = find(abs(H_filt(1:end)) > LIMIT, 1, 'last');
    rangeLow = find(abs(H_filt(1:end)) > LIMIT, 1, 'first');

    if(isempty(rangeHigh))
        rangeHigh = COMP_FLINES;
    end
    if(isempty(rangeLow))
        rangeLow = 1;
    end
    
    fprintf('Hatásos frekvencia sáv: %0.0f - %0.0f Hz \n', f_interp(rangeLow), f_interp(rangeHigh))

    H_mic_temp = H_mic + H_filt;
  
    
    n = 0;
    upperMin = 1;
    lowerMin = 0;
    if (~isempty(find((H_mic_temp(rangeLow:rangeHigh) > tol_interp(1,rangeLow:rangeHigh)) | (H_mic_temp(rangeLow:rangeHigh) < tol_interp(2,rangeLow:rangeHigh)), 1)))
%         lowerMin = 1;
        Ho = H_filt;
        return
    end
    while (abs(upperMin - lowerMin) >= EQUAL_TOLERANCE)
        
        % find minimum difference between upper tolerance band and H_mic_temp
        [upperMin, upperMinPlace] = min(tol_interp(1,rangeLow:rangeHigh)-H_mic_temp(rangeLow:rangeHigh));
        % find minimum difference between lower tolerance band and H_mic_temp
        [lowerMin, lowerMinPlace] = min(H_mic_temp(rangeLow:rangeHigh)-tol_interp(2,rangeLow:rangeHigh));
        
        if ((f_interp(upperMinPlace) > 10100) || (f_interp(lowerMinPlace) > 10100))
            return
        end
        
        if ((lowerMin < 0 || upperMin < 0) && (n == 0))
            return
        end
        
        if((f_interp(lowerMinPlace) > 11500) || (f_interp(upperMinPlace) > 11500))
            return
        end
        
        lowerMin = abs(lowerMin);
        upperMin = abs(upperMin);
        
        dif1 = tol_interp(1,upperMinPlace)- tol_interp(2,upperMinPlace);
        dif2 = tol_interp(1,lowerMinPlace)- tol_interp(2,lowerMinPlace);
        delta_gain = 10;
        
        dir = 1; 
        % if upperMin > lowerMin and gain > 0 => gain has to be bigger
        if ((upperMin > lowerMin) && (gain_dB > 0))
            dir = 1;
        end
        % if upperMin > lowerMin and gain < 0 => abs(gain) has to be bigger
        if ((upperMin > lowerMin) && (gain_dB < 0))
            dir = 1;
        end
        % if upperMin < lowerMin and gain < 0 => abs(gain) has to be smaller
        if ((upperMin < lowerMin) && (gain_dB < 0))
            dir = -1;
        end
        % if upperMin < lowerMin and gain > 0 => gain has to be smaller
        if ((upperMin < lowerMin) && (gain_dB > 0))
            dir = -1;
        end
        gain_dB = gain_dB + dir*delta_gain*1/(2^(n));
        [H_filt,f] = parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
        H_mic_temp = H_mic + H_filt;
        
        %semilogx(f_interp,H_mic_temp)
         % find minimum difference between upper tolerance band and H_mic_temp
        [upperMin, ~] = min(tol_interp(1,rangeLow:rangeHigh)-H_mic_temp(rangeLow:rangeHigh));
        % find minimum difference between lower tolerance band and H_mic_temp
        [lowerMin, ~] = min(H_mic_temp(rangeLow:rangeHigh)-tol_interp(2,rangeLow:rangeHigh));
        
        if (lowerMin < 0 || upperMin < 0)
            gain_dB = gain_dB - dir*delta_gain*1/(2^(n));
        end
        
        lowerMin = abs(lowerMin);
        upperMin = abs(upperMin);
        
        n = n+1;
        if (n == MAX_CYCLE)
            break
        end
         
        
%         pause
    end
    Ho = H_filt;
    fprintf('szûrõ javítva gain: %0.2f \n',gain_dB-initGain)
%     hold off
end