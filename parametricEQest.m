function [gain_dB,fc,bw_oct,Ho,f] = parametricEQest( initGain_dB, initFc, initBw, fs, H_mic, startF, endF, f_interp, H_trgt, tol_interp, configFile)
%parametricEQest - Estimate the parameters of the filter
%
    load(configFile)
    
    

%     e = sum(abs(H_mic));
    [eA,eT] = errorCalc(H_mic,H_trgt,tol_interp,configFile);
    gain_dB = initGain_dB;
    fc = initFc;
    bw_oct = initBw;
    
    % parametric filter with initial values
    [Ho,~]=parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
    % initial error
%     if(sum(abs(H_mic+Ho)) < e)
%         e = sum(abs(H_mic+Ho));
%     end
    [eA_new,eT_new] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
    if(eT_new <= eT)
        if(eA_new <= eA)
            eA = eA_new;
            eT = eT_new;
        end
    end
    switch(EST_MODE)
        case 1
            allError = zeros(2,INNER_LOOP*3*OUTER_LOOP);
            ind = 1;
            for j = 1:OUTER_LOOP
                % fc calibration: from start freq to end freq (uniform (continuous) distribution)
                for i = 1:INNER_LOOP
%                     freq = logspace(log10(startF),log10(endF),40);
%                     param = freq(randi(length(freq)));
                    param = startF * (endF*1.1/startF)^rand;
                    [Ho,~]=parametricEQ(gain_dB,param,bw_oct,fs,f_interp);
        %             if (sum(abs(H_mic+Ho)) < e)
        %                 e = sum(abs(H_mic+Ho));
        %                 fc = param;
        %             end
                    [eA_new,eT_new] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
                    allError(:,ind) = [eA_new;eT_new];
                    ind = ind + 1;
                    if(eT_new <= eT)
                        if(eA_new <= eA)
                            eA = eA_new;
                            eT = eT_new;
                            fc = param;
                        end
                    end
                end
                % gain calibration: form 0.15*MaxGain to 1.15*MaxGain (MaxGain = initGain)
                for i = 1:INNER_LOOP
%                     param = rand + 0.15;
                    param = 12 * rand - 6;
                    [Ho,~]=parametricEQ(initGain_dB+param,fc,bw_oct,fs,f_interp);
        %             if (sum(abs(H_mic+Ho)) < e)
        %                 e = sum(abs(H_mic+Ho));
        %                 gain_dB = initGain*param;
        %             end
                    [eA_new,eT_new] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
                    allError(:,ind) = [eA_new;eT_new];
                    ind = ind + 1;
                     if((eT_new <= eT) && (eA_new <= eA))
%                         if(eA_new < eA)
                            eA = eA_new;
                            eT = eT_new;
                            gain_dB = initGain_dB+param;
%                         end
                    end
                end
                % bw calibration:  form 0.2*initBw to 1.2*initBw (initBw = bandwidth of the error area)
                for i = 1:INNER_LOOP
                    param = rand + 0.2;
                    [Ho,~]=parametricEQ(gain_dB,fc,initBw*param,fs,f_interp);
        %             if (sum(abs(H_mic+Ho)) < e)
        %                 e = sum(abs(H_mic+Ho));
        %                 bw_oct = initBw*param;
        %             end
                    [eA_new,eT_new] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
                    allError(:,ind) = [eA_new;eT_new];
                    ind = ind + 1;
                    if(eT_new <= eT)
                        if(eA_new <= eA)
                            eA = eA_new;
                            eT = eT_new;
                            bw_oct = initBw*param;
                        end
                    end

                end
            end
        case 2
            for j = 1:CYCLES
                paramSel = randi(3);

                switch(paramSel)
                % fc calibration: from start freq to end freq (uniform (continuous) distribution)
                    case 1
%                         freq = logspace(log10(startF),log10(endF),40);
%                         param = freq(randi(length(freq)));
                        param = startF * (endF/startF)^rand;
%                         fprintf('f=%0.0f, ',param);
                        [Ho,~]=parametricEQ(gain_dB,param,bw_oct,fs,f_interp);
            %             if (sum(abs(H_mic+Ho)) < e)
            %                 e = sum(abs(H_mic+Ho));
            %                 fc = param;
            %             end
                        [eA_new,eT_new] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
                        
                        if(eT_new <= eT)
                            if(eA_new <= eA)
                                eA = eA_new;
                                eT = eT_new;
                                fc = param;
                            end
                        end

                % gain calibration: from 0.15*MaxGain to 1.15*MaxGain (MaxGain = initGain)
                % gain calibration: from -6 to +6 dB of initGain
                    case 2
%                         param = 12 * rand - 6;
                        param = rand + 0.15;
%                         fprintf('A=%0.0f, ',param)
                        [Ho,~]=parametricEQ(initGain_dB*param,fc,bw_oct,fs,f_interp);
            %             if (sum(abs(H_mic+Ho)) < e)
            %                 e = sum(abs(H_mic+Ho));
            %                 gain_dB = initGain*param;
            %             end
                        [eA_new,eT_new] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
                        if(eT_new <= eT)
                            if(eA_new <= eA)
                                eA = eA_new;
                                eT = eT_new;
                                gain_dB = initGain_dB+param;
                            end
                        end

                % bw calibration:  form 0.2*initBw to 1.2*initBw (initBw = bandwidth of the error area)
                    case 3
                        param = rand + 0.2;
                        [Ho,~]=parametricEQ(gain_dB,fc,initBw*param,fs,f_interp);
            %             if (sum(abs(H_mic+Ho)) < e)
            %                 e = sum(abs(H_mic+Ho));
            %                 bw_oct = initBw*param;
            %             end
                        [eA_new,eT_new] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile);
                        if(eT_new <= eT)
                            if(eA_new <= eA)
                                eA = eA_new;
                                eT = eT_new;
                                bw_oct = initBw*param;
                            end
                        end

                end
            end
    end
    % final filter
    [Ho,f]=parametricEQ(gain_dB,fc,bw_oct,fs,f_interp);
%     fprintf('Hibaterület: %d;  Toleranciasávot átlépõ pontok: %d \n', eA,eT)
%     figure()
%     plot(allError(1,:),'x')
%     hold on
%     plot(allError(2,:),'o')
%     plot(eA*ones(1,length(allError)))
%     plot(eT*ones(1,length(allError)))
%     vLines = 0:20:length(allError);
%     for j = 2:length(vLines)-1
%         line([vLines(j) vLines(j)], [0 max(allError(1,:))])
% %         [minT, minPlace] = min(allError(2,vLines(j):vLines(j+1)));
% %         plot(minPlace+(j-1)*20,allError(1,minPlace),'o')
% %         plot(minPlace+(j-1)*20,minT,'x')
%     end
%     title('Errors')
%     legend('hibaterület','toleranciahiba','végeleges hibaterület', 'végleges toleranciasávon kívüli pontok száma')
%     ylabel('Hiba')
%     figure()
%     plot(allError(1,:),allError(2,:),'o')
%     hold on
%     h=plot(eA,eT,'x','LineWidth',8);
%     set(h,'linewidth',5)
%     set(h,'markersize',14)
%     legend('lépések hibája','végeleges hiba')
%     xlabel('Hibaterület')
%     ylabel('Toleranciasávon kívüli pontok száma')
%     pause
   
end