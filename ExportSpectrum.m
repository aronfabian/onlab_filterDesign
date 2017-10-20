PLOT = 1;   % legyen-e ábrázolás (=1) vagy csak számolja és mentse a (pl. terc-) sávos spektrumokat
spectr_felbontas = '1/3';   % Ilyen felbontásra számolunk (lehetséges értékei: '1/1', '1/2', '1/3', '1/6', '1/12', '1/24')

measNames = {
    'CutOf_iPhone_5S_MFRecorder_2017-04-27 17-08-07_SweepOnly_h'
    'CutOf_Samsung Galaxy Ace 20170427_173930_SweepOnly_h'
    'CutOf_Samsung Galaxy Alpha - 20170427_183001_SweepOnly_h'
    'CutOf_Samsung Galaxy Mini 2 - 20170427_174623_SweepOnly_h'
    };

for i=1:length(measNames)
    load([measNames{i} ]);  % impulzus válasz és (fázishelyes) átviteli függvény betöltése
    f=fs*(0:(nPoints/2-1))/nPoints;           % frekvenciatengely

    % tercsávossá alakítás (=ez egyfajta simítás), amít így kapunk az egy durva amplitúdómenet 
    [phone_bX,phone_bFreqs]=DA_narrow2band_spectrum(H(1:length(f),:),f,spectr_felbontas,'-frf');
    
    % számolt fügvény lementése
    save([measNames{i} '_bSpectr_'], 'phone_bX', 'phone_bFreqs', 'spectr_felbontas');
    
    % megjelenítés
    if PLOT
        figure(i)
        semilogx(phone_bFreqs, 20*log10(phone_bX),'bx-')
        title(measNames{i}, 'Interpreter', 'none')
        hold on
        % -------------------------------------------------------------------
        % Ezzel tetszõleges felbontásra ki tudod számolni a frekvenciamenetet
        f_interp=100:1:10000;
        X_interp = interp1(phone_bFreqs, phone_bX, f_interp, 'pchip');
        % -------------------------------------------------------------------
        semilogx(f_interp, 20*log10(X_interp), 'r')
        hold off
        pause
    end
end


