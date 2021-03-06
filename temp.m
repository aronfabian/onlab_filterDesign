clear all;
close all;

START_FREQ = 100; % vizsg�lati freki tartom�ny alja
END_FREQ = 10000; % vizsg�lati freki tartom�ny teteje
COMP_FLINES = 1000; % az �tviteli fgv sz�m�t�sakor haszn�lt felbont�s
PLOT_FLINES = 1000; % a megjelen�t�shez haszn�lt felbont�s

A_FILT = 1; % compensation with A-filter
C_FILT = 1; % compensation with C-filter
NO_FILT = 1; % no compensation

HPF_BUTTER = 1; % use highpass butterworth filter 
LPF_BUTTER = 1; % use lowpass butterworth filter

PAUSE = 0; % 0=run without pausing, 1 = run with pausing

fileNames = {
 'CutOf_iPhone_5S_MFRecorder_2017-04-27 17-08-07_SweepOnly_h_bSpectr'
'CutOf_Samsung Galaxy Ace 20170427_173930_SweepOnly_h_bSpectr'
% 'CutOf_Samsung Galaxy Alpha - 20170427_183001_SweepOnly_h_bSpectr'
%'CutOf_Samsung Galaxy Mini 2 - 20170427_174623_SweepOnly_h_bSpectr'
};

main(START_FREQ,END_FREQ,COMP_FLINES,PLOT_FLINES,A_FILT,C_FILT,NO_FILT,HPF_BUTTER,LPF_BUTTER,PAUSE,fileNames);

function main(START_FREQ,END_FREQ,COMP_FLINES,PLOT_FLINES,A_FILT,C_FILT,NO_FILT,HPF_BUTTER,LPF_BUTTER,PAUSE,fileNames)
    

    f_interp = logspace(log10(START_FREQ),log10(END_FREQ),COMP_FLINES);
    f_interp_plot = logspace(log10(START_FREQ),log10(END_FREQ),PLOT_FLINES);
    
    if (A_FILT == 1)
        load('A_weight.mat');
        A_interp = interp1(A_weight(:,1), A_weight(:,2) , f_interp, 'pchip');
    end
    if (C_FILT == 1)
        load('C_weight.mat');
        C_interp = interp1(C_weight(:,1), C_weight(:,2) , f_interp, 'pchip');
    end
    load('tolerance.mat');
    tol_interp = interp1(tolerance(:,1), tolerance(:,2:3) , f_interp, 'pchip');
    tol_interp_plot = interp1(tolerance(:,1), tolerance(:,2:3) , f_interp_plot, 'pchip');
    
    fileID = fopen('output.txt','w');
    fprintf(fileID,'Results:\n');
    if (isempty(fileNames))
        disp('NO INPUT FILE!')
        fprintf(fileID,'NO INPUT FILE!');
        fclose(fileID);
        return
    end

    for fileNum = 1:length(fileNames)
        fprintf(fileID,'-------------------------------------------\n');
         fprintf(fileID,'%s \n',fileNames{fileNum});
        % filterType = 1 --> calc with A_Filter
        % filterType = 2 --> calc with C_Filter
        % filterType = 3 --> calc with NO_Filter
        for filterType = 1:3
            if ((filterType == 1) && (A_FILT == 0))
                filterType = 2;
            end
            if ((filterType == 2) && (C_FILT == 0))
                filterType = 3;
            end
            if ((filterType == 3) && (NO_FILT == 0))
                break;
            end
            
            load([fileNames{fileNum}]);

            X_interp = interp1(phone_bFreqs, phone_bX, f_interp, 'pchip');
            X_interp_plot = interp1(phone_bFreqs, phone_bX, f_interp_plot, 'pchip');

            % H_mic to dB scale
            H_mic = 20*log10(X_interp);
            
            if((filterType == 1) && (A_FILT == 1))
                % H_mic * A_weight
                H_mic = H_mic - A_interp;
                fprintf(fileID,'\tA weighted\n');
            end
            if((filterType == 2) && (C_FILT == 1))
                % H_mic * A_weight
                H_mic = H_mic - C_interp;
                fprintf(fileID,'\tC weighted\n');
            end
            if(filterType == 3)
                fprintf(fileID,'\tWithout weighting filter\n');
            end

            % find 1kHz
            trgt_f = f_interp(1:end-1)<1000 & f_interp(2:end)>1000;
            % H_mic offset with Magnitude[dB] at 1kHz 
            H_mic = H_mic - H_mic(trgt_f) ;


            % H_mic sz�mol�sa megjelen�t�shez
            H_mic_plot = interp1(f_interp, H_mic, f_interp_plot, 'pchip');

            H_trgt = zeros(length(f_interp),1);
            H_trgt_plot = zeros(length(f_interp_plot),1);
            figure()
            semilogx(f_interp_plot,H_mic_plot)
            title(fileNames{fileNum},'Interpreter', 'none')
            xlabel('Frequency (Hz)')
            ylabel('Magnitude (dB)')
            hold on
            semilogx(f_interp_plot,H_trgt_plot)
            semilogx(f_interp_plot,tol_interp_plot(:,1:2),'r--')
            leg = {'H_{mic}original';'H_{trgt}';'tolerance band';'';'Butterworth filters'};
            legend(leg(1:3))

            H_mic_origin = H_mic;

            % Butterworth Filters
            if ((HPF_BUTTER == 1) || (LPF_BUTTER == 1))
                [Hb1,Hb2,fcb1,fcb2] = butterworthFilters(H_mic,f_interp,START_FREQ,END_FREQ);
                if(HPF_BUTTER == 0)
                    Hb1 = zeros(1,length(f_interp));
                else
                    fprintf(fileID,'\t\tButterworth HPF cutoff freq: %0.1f Hz\n',fcb1);
                end
                if(LPF_BUTTER == 0)
                    Hb2 = zeros(1,length(f_interp));
                else
                    fprintf(fileID,'\t\tButterworth LPF cutoff freq: %0.1f Hz\n',fcb2);
                end
                H_mic = H_mic+Hb1+Hb2;
                
            end
            
            H_mic_plot = interp1(f_interp, H_mic, f_interp_plot, 'pchip');
            semilogx(f_interp_plot,H_mic_plot)
            legend(leg(1:5))

            if(PAUSE == 1)
                pause
            end
            
            filterNum = 0;
            finishFlag = 0;
            while(finishFlag == 0)
                filterNum = filterNum + 1;
                finishFlag = 0;
                
                % fine-tuning of prev. filters
                if (filterNum > 2)
                    for i = 1:filterNum-1
                        [H1,~] = parametricEQ(filters(i,1),filters(i,2),filters(i,3),fs,f_interp);
                        [filters(i,1),filters(i,2),filters(i,3),H2,f] = parametricEQest(filters(i,1),filters(i,2),filters(i,3),fs,H_mic-H1,filters(i,4),filters(i,5),f_interp);
                        H_mic = H_mic - H1 + H2;
                    end
                end

                % H_mic and H_trgt=0 crossings
                cross = find((H_mic(1:end-1) .* H_mic(2:end) < 0));
                areaNum = length(cross)+1;

                % start-end frequency
                startEndFreq = [START_FREQ f_interp(cross) END_FREQ];
                startEnd = [1 cross length(f_interp)];

                % error areas
                errorAreas = zeros(areaNum,1);
                error = H_mic' - H_trgt;
                for i = 1:areaNum
                    errorAreas(i) = sum(abs(error(startEnd(i):startEnd(i+1))));
                end

                % maximum errors
                maxErrors = zeros(areaNum,1);
                maxErrorsFreq = zeros(areaNum,1);
                maxErrorsPlace = zeros(areaNum,1);
                for i = 1:areaNum
                    if (H_mic(startEnd(i)+1) > 0)
                        [maxErrors(i), maxErrorsPlace(i)] = max(H_mic(startEnd(i):startEnd(i+1)));
                    else
                        [maxErrors(i), maxErrorsPlace(i)] = min(H_mic(startEnd(i):startEnd(i+1)));
                    end
                    maxErrorsPlace(i) = maxErrorsPlace(i) + startEnd(i)-1;
                    maxErrorsFreq(i) =  f_interp(maxErrorsPlace(i));
                    % if a maximum error of an error area is inside the tolerance
                    % band, error area will be 0 (it's good enough)
                    if (isempty(find((H_mic(startEnd(i):startEnd(i+1)) > tol_interp(startEnd(i):startEnd(i+1),1)) | (H_mic(startEnd(i):startEnd(i+1)) < tol_interp(startEnd(i):startEnd(i+1),2)), 1)))
                         errorAreas(i) = 0;
                    end
                end

                % end of iteration
                if(max(errorAreas) == 0)
                    fprintf('Class1-es tolerancia s�vba tartoz� �tvitelhez sz�ks�ges \nsz�r�k sz�ma: ')
                    disp(filterNum-1)
                    disp('-------------------------')
                    figure
                    H_mic_o_plot = interp1(f_interp, H_mic_origin, f_interp_plot, 'pchip');
                    semilogx(f_interp_plot,H_mic_o_plot)
                    title('Original and final transfer, Filters')
                    hold on
                    semilogx(f_interp_plot,tol_interp_plot(:,1:2),'r--')
                    legend( [leg(1); leg(3:4)])
                    Hb1_plot = interp1(f_interp,Hb1,f_interp_plot,'pchip');
                    Hb2_plot = interp1(f_interp,Hb2,f_interp_plot,'pchip');
                    semilogx(f_interp_plot,Hb1_plot,'DisplayName', 'Butterworth HPF')
                    semilogx(f_interp_plot,Hb2_plot,'DisplayName', 'Butterworth LPF')
                    H_mic_o_plot = H_mic_o_plot+Hb1_plot+Hb2_plot;


                    for i = 1:size(filters,1)
                        [He,~] = parametricEQ(filters(i,1),filters(i,2),filters(i,3),fs,f_interp_plot);
                        H_mic_o_plot = H_mic_o_plot + He; 
                        semilogx(f_interp_plot,He,'DisplayName', string(i))
                        
                        fprintf(fileID,'\t\tParametric Filter %d: gain: %0.2f dB, fc: %0.1f Hz, bandwidth: %0.4f octave\n',i,filters(i,1),filters(i,2),filters(i,3));

                    end

                    semilogx(f_interp_plot,H_mic_o_plot,'DisplayName', 'H_{final}')       
                    finishFlag = 1;

                end
                
                if(finishFlag == 0)
                    % select maximum error area
                    [~, maxAreaNum] = max(errorAreas);
                    % initial center frequnecy: frequency of the maxumim magnitude in
                    % the selected area
                    fc=maxErrorsFreq(maxAreaNum);
                    % initial bandwidth in octave: bandwidth of the selected area
                    bw=log2(startEndFreq(maxAreaNum+1)/startEndFreq(maxAreaNum));
                    % initial gain: (-1)*maximum magnitude in the selected area
                    gain = (-1)*maxErrors(maxAreaNum);
                    fs=44100;

                    fprintf('Frekvencia s�v: %0.0f - %0.0f Hz \n', startEndFreq(maxAreaNum), startEndFreq(maxAreaNum+1))
                    fprintf('fc: %0.0f Hz, bw: %0.2f oct, gain: %0.2f dB \n', fc, bw, gain)

                    %pause
                    % estimate parametric filter
                    [estGain,estFc,estBw,Ho,~] = parametricEQest(gain,fc,bw,fs,H_mic,startEndFreq(maxAreaNum),startEndFreq(maxAreaNum+1),f_interp);
                    [estGain,estFc,estBw,Ho,~] = toleranceOptim(estGain,estFc,estBw,fs,H_mic,f_interp,COMP_FLINES, tol_interp);
                    [Ho,f] = parametricEQ(estGain,estFc,estBw,fs,f_interp);
                    semilogx(f,Ho,'DisplayName', string(filterNum))
                    % save filter parameters 
                    if (filterNum == 1)
                        filters =  [estGain,estFc,estBw,startEndFreq(maxAreaNum),startEndFreq(maxAreaNum+1)];
                    else
                        filters = [filters; estGain,estFc,estBw,startEndFreq(maxAreaNum),startEndFreq(maxAreaNum+1)];
                    end

                    disp('A becsl� �rt�kei: ')
                    fprintf('estFc: %0.0f Hz, estBw: %0.2f oct, estGain: %0.2f dB \n', estFc , estBw, estGain)
                    disp('-------------------------')


                    % new transfer function = transfer function + parametric filter
                    H_mic = H_mic + Ho;
                    H_mic_plot = interp1(f_interp, H_mic, f_interp_plot, 'pchip');

                    % plot new transfer function
                    semilogx(f_interp_plot, H_mic_plot, 'DisplayName', string(filterNum))

                    if(PAUSE == 1)
                        pause
                    end
                end
            end
        end
    end % fileNum
    
    fclose(fileID);
    
end
