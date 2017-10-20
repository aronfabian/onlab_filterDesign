function varargout =DA_narrow2band_spectrum(xSpectr,nFreqs,OctOrThird,varargin)
% [bX,bFreqs]=da_narrow2band_spectrum(xSpectr,nFreqs,OctOrThird,varargin)
%
% Data Acoustics - Converter routine
%
% Generates Octave or Third octave, or finer band spectra from narrow band spectra.
% <OctOrThird> must be one of the followings: 
%    '1/1', '1/2', '1/3', '1/6', '1/12', '1/24'
%
%
% xSpectr:      Matrix, with spectras of the channels in its COLUMNs
% nFreqs:       Column, with the frequency values of the spectral lines
% OctOrThird:   must be one of the followings: '1/1', '1/2', '1/3', '1/6', '1/12', '1/24'
%
% varargin:
%   '-frf':     specify this, if your input is a transfer function, rahter then an autospectrum
%               (In this case each band's value must be divided by the number of spectral lines
%               contributing to the given band's value. This option lets perform the correction.)
%   '-lpb' <n>: Spectral band is valid if contains at least <n> spectral lines. Default = 5.
%
% bX:           Matrix, contains the Octave, Third Octave, etc. band spectras generated from xSpectr.
%               (linear values, not dBs)
% bFreqs:       Frequency values for each line of the spectras.
%
% Def: X(band) = sqrt( szum( X(k) * conj(X(k)) ) )
% A.B. Nagy, F. Marki - 2003-2011.


% 'conj' beszúrása az 51-es sorba

prg_ver='2.1';

ERR = 'DA_narrow2band_spectrum error';

if ischar(xSpectr)
  if  strcmp(xSpectr,'-rms_spectr')
      if strcmp(nFreqs(end-10:end),'_rms_spectr')
          matfile_spectr=[nFreqs(1:end-11) '_3rd_spectr.mat'];
      else
          matfile_spectr=[nFreqs '_3rd_spectr.mat'];
      end
      matfile = [nFreqs '.mat'];
  else
      matfile = [xSpectr '_temp.mat'];
  end
  clear xSpectr;
  
  try
    if ~exist(matfile,'file')                                               % check the file
        error(['Frequency data file (' matfile ') not be found!']);
    end
 %   hWaitbar=waitbar(0,'Reading frequency data file...');
%    set(hWaitbar,'Name',['DA_narrow2band_spectrum - ' prg_ver]);
    load('-mat', matfile);                                                  % load database
    if ~exist('X','var')
%        close(hWaitbar);
        matfile
        error('Invalid file specified. Frequency data doesn''t not exist!');
    end
    xSpectr=X;
  catch
    errordlg(lasterr, ERR);        
  end
else
    if isstruct(nFreqs)
        info=nFreqs;
        nFreqs=info.nFreqs;
    end
 %   hWaitbar=waitbar(0,'Computing narrow band spectra...');
%    set(hWaitbar,'Name',['DA_wav2spect - ' prg_ver]);
end    


bX=[]; bFreqs=[]; start_fline=[];

isFRF=0;
energyCORR=1;
min_lines_per_band=5;

isDataToProcess=1;
for n=1:nargin-3
    if isDataToProcess
        switch lower(varargin{n}) 
            case '-frf', isFRF=1;
                         disp('FRF correction will be applied.');
            case '-energy_corr', if nargin-3 < n+1
                                     disp(['ERROR: Energy correction value doesn''t exist']); 
                                     return;
                                 end
                                 energyCORR=varargin{n+1};
                                 if ~isnumeric(energyCORR)
                                     disp(['ERROR: Energy correction should be a numeric value']); 
                                     return;
                                 end
                                 isDataToProcess=0;
            case '-lpb', if nargin-3 < n+1
                                     disp(['ERROR: Lines Per Band data doesn''t exist']); 
                                     return;
                                 end
                                 min_lines_per_band=varargin{n+1};
                                 if ~isnumeric(min_lines_per_band) | round(min_lines_per_band)~=min_lines_per_band | min_lines_per_band<=0
                                     disp(['ERROR: Lines Per Band data should be positive integer']); 
                                     return;
                                 end
                                 isDataToProcess=0;
                         
            otherwise disp(['WARNING: Unrecognized option ''' varargin{n} '''']);
        end
    else
        isDataToProcess=1;
    end
end

if ~exist('isFRF','var'),   isFRF=0;  end

if min(size(nFreqs))~=1
    disp(['ERROR: Frequency data must be a row or coloumn VECTOR!']); 
    return;
end
if size(xSpectr,1)~=length(nFreqs)
    xSpectr=xSpectr.';
    disp(['WARNING: Spectrum data should be a coloumn vector or coloumns of a matrix']); 
end
if size(xSpectr,1)~=length(nFreqs)
    disp(['ERROR: Length of frequency vector must be equal to the size of spectrum/-a vector(s)']); 
    return;
end

if nFreqs(1)==0 xSpectr(1,:)=[]; nFreqs(1)=[]; end

tempCFreqs12=[16 22.4 31.5 45 63 90 125 180 250 355 500 710 1000 1400 2000 2800 4000 5600 8000 11200 16000];
tempCFreqs13=[1 1.25 1.6 2 2.5 3.15 4 5 6.3 8 10];

n_divide=-1;
tCF=tempCFreqs13;
switch OctOrThird
    case '1/1',  tCF=tempCFreqs12(1:2:end);
    case '1/2',  tCF=tempCFreqs12;
    case '1/3',  n_divide=0;
    case '1/6',  n_divide=1;
    case '1/12', n_divide=2;
    case '1/24', n_divide=3;
    otherwise disp(['ERROR: unknown bandwidth: ''' OctOrThird '''']);   return;
end
for n=1:n_divide
    temp= sqrt(tCF(1:end-1).*tCF(2:end));
    tCF(1:2:end*2-1)=tCF;
    tCF(2:2:end-1)=temp;
end
tCF=tCF(1:end-1);

min_kitevo=floor(log10(min(nFreqs))/3-1)*3; % mivel az oktávsáv 3 dekádonként ismétlõdik
max_kitevo=ceil(log10(max(nFreqs)));        % ha pont kerek dekád a max. frek. akkor is egy dekáddal fölé megyünk

tempCFreqs=[];
for n=min_kitevo:max_kitevo
    if n_divide>=0 | mod(n,3)==0
        tempCFreqs=[tempCFreqs  tCF*10^n];
    end
end
% ---------------------- tercsáv középfrekvenciák elkészültek ---------------------------

[nFreqLines,nCh]=size(xSpectr);
bX=zeros(length(tempCFreqs),nCh);
only4rms=zeros(length(tempCFreqs),1);

xLimits=sqrt(tempCFreqs(1:end-1).*tempCFreqs(2:end));

band_no=1; start_fline=0;
for f=1:nFreqLines,
%     if nFreqs(f) > xLimits(band_no)
        while nFreqs(f) > xLimits(band_no),
%             if start_fline(band_no)<f-1 && start_fline(band_no)>0
% [                xSpectr(start_fline(band_no):f-1,:) ; only4rms(band_no)]
%                 pause
%             end
            band_no=band_no+1;
            start_fline(band_no)=f;
        end
%     end
%    band_no
%    pause
    % A következõ sorban kell a conj, hogy komplex spektrumokat is lehessen használni    
    % Negyzetosszeg generalasa egy savra, az összes spektrumra külön-külön
    bX(band_no,:)=bX(band_no,:)+(xSpectr(f,:).*conj(xSpectr(f,:)));
    % Nem kene itten osztani a savba besorolt mintak szamaval, hogy tenyleg
    % rms legyen???
    % NEM, csak FRF-ek esetén, azt pedig késõbb vesszük figyelembe

    % A kovetkezoben szamolom, hogy egy savba hany spektrumvonal esett
    only4rms(band_no)=only4rms(band_no)+1;
end

%[(1:length(tempCFreqs))' tempCFreqs' only4rms only4rms<min_lines_per_band]

%min_lines_per_band
%only4rms
valid_pos=find(only4rms>=min_lines_per_band);    % Megkeressük az érvényes frekvenciasávokat
last_valid_band=max(valid_pos);
only4rms(last_valid_band+1:end)=[];
not_valid_pos=find(only4rms<min_lines_per_band); % Mostmár csak az érvényes tartomány alatt vannak érvénytelen sávok
first_valid_band=max(not_valid_pos)+1;
only4rms(1:first_valid_band-1)=[];

start_fline=start_fline(first_valid_band);

bX=bX(first_valid_band:last_valid_band,:);
bFreqs=tempCFreqs(first_valid_band:last_valid_band)';

% [bFreqs bX only4rms]
% Ha FRF-rõl van szó, akkor le kell osztani a sávba esõ spektrumvonalak számával
if isFRF
    for k=1:length(only4rms)
        bX(k,:)=bX(k,:)/only4rms(k);    
    end
end
% Es vegul johet a gyokvonas
if energyCORR~=1
    disp(['Energy correction will be applied: ' num2str(energyCORR) ]);
    bX=sqrt(bX*energyCORR);
else
    bX=sqrt(bX);
end
% [bFreqs bX only4rms]

if exist('info','var')
    info.History(end+1) = { ['DA_narrow2band_spectrum. Bandwidth: ' OctOrThird]  };
end
if nargout>1
    varargout{1}=bX;
    if exist('info','var')
        info.bFreqs=bFreqs;
        varargout{2}=bFreqs;
        varargout{3}=info;
    else
        varargout{2}=bFreqs;
        varargout{3}=start_fline;
    end
    
else
%    info.bandfreq=freq;
%    waitbar(1, hWaitbar, 'Saving results...');
    save(matfile_spectr,'info','bFreqs','bX');
    varargout{1}=info;
end
%close(hWaitbar);
