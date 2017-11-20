function [gain_dB,fc,bw_oct,Ho,f] = toleranceOptim( initGain, initFc, initBw, fs, H_mic, f_interp,COMP_FLINES, tol_interp)
%toleranceOptim - modifify the parameters of the parametric filter
%therefore the tansfer function is in the middle of the tolerance band
%

    gain_dB = initGain;
    fc = initFc;
    bw_oct = initBw;

    % effective range of the filter
    limit = 0.2;
    [H_filt,f] = parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
    rangeHigh = find(H_filt(1:end-1)>limit & H_filt(2:end)<limit);
    rangeLow = find(H_filt(1:end-1)<limit & H_filt(2:end)>limit);
    if(isempty(rangeHigh))
        rangeHigh = COMP_FLINES;
    end
    if(isempty(rangeLow))
        rangeLow = 1;
    end
    
    H_mic_temp = H_mic + H_filt;
    
    % goal: lowerMin = upperMin
    
    % changing of transfer function at upperMinPlace and lowerMinPlace =>
    % to calcultate new gain
    
    % OR ~ successive approximation
%     figure
%     semilogx(f_interp,H_mic_temp)
%     hold on
%     semilogx(f_interp,tol_interp(:,1:2),'r--')
    
    n = 0;
    upperMin = 1;
    lowerMin = 0;
    if (~isempty(find((H_mic_temp(rangeLow:rangeHigh) > tol_interp(rangeLow:rangeHigh,1)) | (H_mic_temp(rangeLow:rangeHigh) < tol_interp(rangeLow:rangeHigh,2)), 1)))
%         lowerMin = 1;
        Ho = H_filt;
        return
    end
    while (abs(upperMin - lowerMin) >= 0.2)
        
        % find minimum difference between upper tolerance band and H_mic_temp
        [upperMin, upperMinPlace] = min(tol_interp(rangeLow:rangeHigh,1)-H_mic_temp(rangeLow:rangeHigh)');
        % find minimum difference between lower tolerance band and H_mic_temp
        [lowerMin, lowerMinPlace] = min(H_mic_temp(rangeLow:rangeHigh)'-tol_interp(rangeLow:rangeHigh,2));
        
        if ((lowerMin < 0 || upperMin < 0) && (n == 0))
            return
        end
        
        lowerMin = abs(lowerMin);
        upperMin = abs(upperMin);
        
        dif1 = tol_interp(upperMinPlace,1)- tol_interp(upperMinPlace,2);
        dif2 = tol_interp(lowerMinPlace,1)- tol_interp(lowerMinPlace,2);
        if(dif1 > dif2)
            delta_gain = 10;
        else
            delta_gain = 10;
        end
        
        dir = 1; % dir = 1 => bigger, dir = -1 => smaller
        % if upperMin > lowerMin and gain > 0 => gain has to be bigger
        if ((upperMin > lowerMin) && (gain_dB > 0))
            dir = 1;
        end
        % if upperMin > lowerMin and gain < 0 => gain has to be smaller
        if ((upperMin > lowerMin) && (gain_dB < 0))
            dir = 1;
        end
        % if upperMin < lowerMin and gain < 0 => gain has to be bigger
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
        [upperMin, ~] = min(tol_interp(rangeLow:rangeHigh,1)-H_mic_temp(rangeLow:rangeHigh)');
        % find minimum difference between lower tolerance band and H_mic_temp
        [lowerMin, ~] = min(H_mic_temp(rangeLow:rangeHigh)'-tol_interp(rangeLow:rangeHigh,2));
        
        if (lowerMin < 0 || upperMin < 0)
            gain_dB = gain_dB - dir*delta_gain*1/(2^(n));
        end
        
        lowerMin = abs(lowerMin);
        upperMin = abs(upperMin);
        
        n = n+1;
        if (n == 1000)
            break
        end
         
        
%         pause
    end
    Ho = H_filt;
    fprintf('szûrõ javítva gain: %0.2f \n',gain_dB-initGain)
%     hold off
end