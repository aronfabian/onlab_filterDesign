% close all
%clear all

%config();
main('config.mat');

function main(configFile)

    load(configFile)
    
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
        
        if (TOLERANCE == 1)
            from = 1;
            to = 1;
        end
        if (TOLERANCE == 2)
            from = 2;
            to = 2;
        end
        if (TOLERANCE == 3)
            from = 1;
            to = 2;
        end
        for tol_index = from:to
            fprintf(fileID,'Class%d:\n', tol_index);
            tol_interp = (interp1(tolerance(:,1), tolerance(:,tol_index*2:tol_index*2+1) , f_interp, 'pchip'))';
            tol_interp_plot = (interp1(tolerance(:,1), tolerance(:,tol_index*2:tol_index*2+1) , f_interp_plot, 'pchip'))';
            
            % filterType = 1 --> calc with A_Filter
            % filterType = 2 --> calc with C_Filter
            % filterType = 3 --> calc with NO_Filter
            filterType = 1;
            while (filterType < 4)
                if ((filterType == 1) && (A_FILT == 0))
                    filterType = 2;
                end
                if ((filterType == 2) && (C_FILT == 0))
                    filterType = 3;
                end
                if ((filterType == 3) && (NO_FILT == 0))
                    break;
                end

                figTitle = 0;
                switch filterType
                    case 1
                        fprintf('A szûrõ\n');
                        figTitle = 'A weighted';
                    case 2
                        fprintf('C szûrõ\n');
                        figTitle = 'C weighted';
                    case 3
                        fprintf('nincs szûrõ\n');
                        figTitle = 'No weighting filter';
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
                    % H_mic * C_weight
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


                % H_mic számolása megjelenítéshez
                H_mic_plot = interp1(f_interp, H_mic, f_interp_plot, 'pchip');

                H_trgt = zeros(1,length(f_interp));
                H_trgt_plot = zeros(1,length(f_interp_plot));
                figure()
                semilogx(f_interp_plot,H_mic_plot)
                title({fileNames{fileNum};['Class' num2str(tol_index)]; figTitle},'Interpreter', 'none')
                xlabel('Frequency (Hz)')
                ylabel('Magnitude (dB)')
                hold on
                semilogx(f_interp_plot,H_trgt_plot)
                semilogx(f_interp_plot,tol_interp_plot(1:2,:),'r--')
                leg = {'H_{mic}original';'H_{trgt}';'tolerance band';''};
                legend(leg(1:3))

                H_mic_origin = H_mic;
                legend(leg(1:4))

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
                            if(filters(i,6) == 0)
                                [H1,~] = parametricEQ(filters(i,1),filters(i,2),filters(i,3),fs,f_interp);
                                [filters(i,1),filters(i,2),filters(i,3),H2,f] = parametricEQest(filters(i,1),filters(i,2),filters(i,3),fs,H_mic-H1,filters(i,4),filters(i,5),f_interp,H_trgt,tol_interp,configFile);
                                H_mic = H_mic - H1 + H2;
                                fprintf('Finomhangolás (%d): estFc: %0.0f Hz, estBw: %0.2f oct, estGain: %0.2f dB \n', i, filters(i,2) , filters(i,3), filters(i,1))

                            end
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
    %                 error = H_mic - H_trgt;
                    for i = 1:areaNum
    %                     errorAreas(i) = sum(abs(error(startEnd(i):startEnd(i+1))));
                        [errorAreas(i),~] = errorCalc(H_mic(startEnd(i):startEnd(i+1)),H_trgt(startEnd(i):startEnd(i+1)),tol_interp(:,startEnd(i):startEnd(i+1)),configFile);
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
                        if (isempty(find((H_mic(startEnd(i):startEnd(i+1)) > tol_interp(1,startEnd(i):startEnd(i+1))) | (H_mic(startEnd(i):startEnd(i+1)) < tol_interp(2,startEnd(i):startEnd(i+1))), 1)))
                             errorAreas(i) = 0;
                        end
                    end

                    % end of iteration
                    if(max(errorAreas) == 0)
                        fprintf('Class%d-es tolerancia sávba tartozó átvitelhez szükséges \nszûrõk száma: ',tol_index)
                        disp(filterNum-1)
                        disp('-------------------------------------------------')
                        figure
                        H_mic_o_plot = interp1(f_interp, H_mic_origin, f_interp_plot, 'pchip');
                        semilogx(f_interp_plot,H_mic_o_plot)
                        title({'Original and final transfer, Filters';['Class' num2str(tol_index)];figTitle})
                        hold on
                        semilogx(f_interp_plot,tol_interp_plot(1:2,:),'r--')
                        legend( [leg(1); leg(3:4)])
                           
                        for i = 1:size(filters,1)
                            if (filters(i,6) == 0) % parametric filter
                                [He,~] = parametricEQ(filters(i,1),filters(i,2),filters(i,3),fs,f_interp_plot);
                                H_mic_o_plot = H_mic_o_plot + He; 
                                semilogx(f_interp_plot,He,'DisplayName', string(i))

                                fprintf(fileID,'\t\tParametric Filter %d: gain: %0.2f dB, fc: %0.1f Hz, bandwidth: %0.4f octave\n',i,filters(i,1),filters(i,2),filters(i,3));

                            else % butterworth filter
                                if(filters(i,6) == 1)
                                    [b,a] = butter(1,filters(i,2)/(fs/2),'high');
                                    fprintf(fileID,'\t\tButterworth HPF cutoff freq: %0.1f Hz\n',filters(i,2));
                                else
                                    [b,a] = butter(1,filters(i,2)/(fs/2),'low');
                                    fprintf(fileID,'\t\tButterworth LPF cutoff freq: %0.1f Hz\n',filters(i,2));
                                end
                                [He,~] = freqz(b,a,f_interp_plot,fs);
                                He = 20*log10(abs(He));
                                H_mic_o_plot = H_mic_o_plot + He; 
                                semilogx(f_interp_plot,He,'DisplayName', string(i))
                            end
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

                        fprintf('Frekvencia sáv: %0.0f - %0.0f Hz \n', startEndFreq(maxAreaNum), startEndFreq(maxAreaNum+1))
                        fprintf('Kezdeti értékek: fc: %0.0f Hz, bw: %0.2f oct, gain: %0.2f dB \n', fc, bw, gain)

                        %pause
                        % estimate parametric filter
                        [estGain,estFc,estBw,Ho,~] = parametricEQest(gain,fc,bw,fs,H_mic,startEndFreq(maxAreaNum),startEndFreq(maxAreaNum+1),f_interp,H_trgt,tol_interp,configFile);
                        fprintf('A becslõ értékei: estFc: %0.0f Hz, estBw: %0.2f oct, estGain: %0.2f dB \n', estFc , estBw, estGain)
                        [estGain,estFc,estBw,Ho,~] = toleranceOptim(estGain,estFc,estBw,fs,H_mic,f_interp,COMP_FLINES, tol_interp,configFile);
    %                     [Ho,f] = parametricEQ(estGain,estFc,estBw,fs,f_interp);
    %                     semilogx(f,Ho,'DisplayName', string(filterNum))
                         butterType = 0;
                        %try butterworth filters
                         if ((HPF_BUTTER == 1) || (LPF_BUTTER == 1))
                            [Hb1,Hb2,fcb1,fcb2] = butterworthFilters(H_mic,f_interp,tol_interp,H_trgt,configFile);
                            if(HPF_BUTTER == 0)
                                Hb1 = zeros(1,length(f_interp));
                            else
                                %fprintf(fileID,'\t\tButterworth HPF cutoff freq: %0.1f Hz\n',fcb1);
                            end
                            if(LPF_BUTTER == 0)
                                Hb2 = zeros(1,length(f_interp));
                            else
                                %fprintf(fileID,'\t\tButterworth LPF cutoff freq: %0.1f Hz\n',fcb2);
                            end
                            [hpfErrorA,hpfErrorT] = errorCalc(H_mic+Hb1,H_trgt,tol_interp,configFile);
                            [lpfErrorA,lpfErrorT] = errorCalc(H_mic+Hb2,H_trgt,tol_interp,configFile);
                            [paramErrorA,paramErrorT] = errorCalc(H_mic+Ho,H_trgt,tol_interp,configFile); 
                            [errorA,errorT] = errorCalc(H_mic,H_trgt,tol_interp,configFile); 
                            E = [hpfErrorT lpfErrorT paramErrorT errorT; hpfErrorA lpfErrorA paramErrorA errorA ];
                            [minT,minPlace] = min(E(1,:));
                            for n = 1:4
                                if ((minT == E(1,n)) && E(2,minPlace) > E(2,n))
                                    minPlace = n;
                                    minT = E(2,n);
                                end
                            end
                            switch (minPlace)
                                case 1
                                    estGain = 0;
                                    estFc = fcb1;
                                    estBw = 0;
                                    butterType = 1; %HPF
                                case 2 
                                    estGain = 0;
                                    estFc = fcb2;
                                    estBw = 0;
                                    butterType = 2; % LPF
                                case 3 
                                    butterType = 0;
                                case 4
                                    butterType = 0;
                                    disp('\n Egyik szûrõvel sem sikerült jobb átvitelt kialakítani!')
                            end


                        end

                        % save filter parameters 
                        if (filterNum == 1)
                            filters =  [estGain,estFc,estBw,startEndFreq(maxAreaNum),startEndFreq(maxAreaNum+1),butterType];
                        else
                            filters = [filters; estGain,estFc,estBw,startEndFreq(maxAreaNum),startEndFreq(maxAreaNum+1),butterType];
                        end

    %                     disp('A becslõ értékei: ')
                        fprintf('Tol. sávba optimalizált: estFc: %0.0f Hz, estBw: %0.2f oct, estGain: %0.2f dB \n', estFc , estBw, estGain)
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
                filterType = filterType + 1;
            end % filterType
        end % tol_index
    end % fileNum

    fclose(fileID);

end

