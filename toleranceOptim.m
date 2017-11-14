function [gain_dB,fc,bw_oct,Ho,f] = toleranceOptim( initGain, initFc, initBw, fs, H_mic, f_interp,COMP_FLINES, tol_interp)
%toleranceOptim - modifify the parameters of the parametric filter
%therefore the tansfer function is in the middle of the tolerance band
%

    gain_dB = initGain;
    fc = initFc;
    bw_oct = initBw;

    % effective range of the filter
    limit = 0.1;
    [H_filt,f] = parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
    rangeHigh = find(H_filt(1:end-1)>limit & H_filt(2:end)<limit);
    rangeLow = find(H_filt(1:end-1)<limit & H_filt(2:end)>limit);
    if(isempty(rangeHigh))
        rangeHigh = COMP_FLINES;
    end
    if(isempty(rangeLow))
        rangeLow = 1;
    end
    range =  [rangeLow rangeHigh];
    
   
    
    % goal: lowerMin = upperMin
    
    % changing of transfer function at upperMinPlace and lowerMinPlace =>
    % to calcultate new gain
    
    % OR ~ successive approximation
    figure
    semilogx(f_interp,H_mic)
    hold on
    semilogx(f_interp,tol_interp(:,1:2),'r--')
    
    n = 0;
    upperMin = 1;
    lowerMin = 0;
    if (~isempty(find((H_mic(rangeLow:rangeHigh) > tol_interp(rangeLow:rangeHigh,1)) | (H_mic(rangeLow:rangeHigh) < tol_interp(rangeLow:rangeHigh,2)))))
        lowerMin = 1;
    end
    while (abs(upperMin - lowerMin) >= 0.2)
        
        % find minimum difference between upper tolerance band and H_mic
        [upperMin, upperMinPlace] = min(tol_interp(rangeLow:rangeHigh,1)-H_mic(rangeLow:rangeHigh)');
        % find minimum difference between lower tolerance band and H_mic
        [lowerMin, lowerMinPlace] = min(H_mic(rangeLow:rangeHigh)'-tol_interp(rangeLow:rangeHigh,2));
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
        if ((upperMin > lowerMin) && (gain_dB > 0))
            dir = -1;
        end
        % if upperMin < lowerMin and gain < 0 => gain has to be bigger
        if ((upperMin > lowerMin) && (gain_dB > 0))
            dir = 1;
        end
        % if upperMin < lowerMin and gain > 0 => gain has to be smaller
        if ((upperMin > lowerMin) && (gain_dB > 0))
            dir = -1;
        end
        gain_dB = gain_dB + dir*delta_gain*1/(2^(n));
        H_mic = H_mic - H_filt;
        [H_filt,f] = parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
        H_mic = H_mic + H_filt;
        n = n+1;
        semilogx(f_interp,H_mic)
         % find minimum difference between upper tolerance band and H_mic
        [upperMin, upperMinPlace] = min(tol_interp(rangeLow:rangeHigh,1)-H_mic(rangeLow:rangeHigh)');
        % find minimum difference between lower tolerance band and H_mic
        [lowerMin, lowerMinPlace] = min(H_mic(rangeLow:rangeHigh)'-tol_interp(rangeLow:rangeHigh,2));
        lowerMin = abs(lowerMin);
        upperMin = abs(upperMin);
        
        pause
    end
    H_mic = H_mic - H_filt;
    Ho = H_filt;
    hold off
end