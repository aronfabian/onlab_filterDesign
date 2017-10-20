PLOT = 1;   % legyen-e �br�zol�s (=1) vagy csak sz�molja �s mentse a (pl. terc-) s�vos spektrumokat
spectr_felbontas = '1/3';   % Ilyen felbont�sra sz�molunk (lehets�ges �rt�kei: '1/1', '1/2', '1/3', '1/6', '1/12', '1/24')

measNames = {
    'CutOf_iPhone_5S_MFRecorder_2017-04-27 17-08-07_SweepOnly_h'
    'CutOf_Samsung Galaxy Ace 20170427_173930_SweepOnly_h'
    'CutOf_Samsung Galaxy Alpha - 20170427_183001_SweepOnly_h'
    'CutOf_Samsung Galaxy Mini 2 - 20170427_174623_SweepOnly_h'
    };

for i=1:length(measNames)
    load([measNames{i} ]);  % impulzus v�lasz �s (f�zishelyes) �tviteli f�ggv�ny bet�lt�se
    f=fs*(0:(nPoints/2-1))/nPoints;           % frekvenciatengely

    % tercs�voss� alak�t�s (=ez egyfajta sim�t�s), am�t �gy kapunk az egy durva amplit�d�menet 
    [phone_bX,phone_bFreqs]=DA_narrow2band_spectrum(H(1:length(f),:),f,spectr_felbontas,'-frf');
    
    % sz�molt f�gv�ny lement�se
    save([measNames{i} '_bSpectr_'], 'phone_bX', 'phone_bFreqs', 'spectr_felbontas');
    
    % megjelen�t�s
    if PLOT
        figure(i)
        semilogx(phone_bFreqs, 20*log10(phone_bX),'bx-')
        title(measNames{i}, 'Interpreter', 'none')
        hold on
        % -------------------------------------------------------------------
        % Ezzel tetsz�leges felbont�sra ki tudod sz�molni a frekvenciamenetet
        f_interp=100:1:10000;
        X_interp = interp1(phone_bFreqs, phone_bX, f_interp, 'pchip');
        % -------------------------------------------------------------------
        semilogx(f_interp, 20*log10(X_interp), 'r')
        hold off
        pause
    end
end


