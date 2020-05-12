% A-weight approximator

clear
fs=44100;
db = jsondecode(fileread('PhoneDatabase__pretty.json'));
uncalDev = db.Phones(find(strcmp({db.Phones.phone_id},'uncalibrated_device')));

%%
load A_weight

EQs = uncalDev.Class1.A_weighted;

H=[];
for i=1:size(EQs.HPF)
    [H(:,end+1),f]=ResonantFilter(EQs.HPF(i).fc,'HPF',fs);
end
for i=1:size(EQs.LPF)
    [H(:,end+1),f]=ResonantFilter(EQs.LPF(i).fc,'LPF',fs);
end
for i=1:size(EQs.Parametric)
    filt=EQs.Parametric(i);
    [H(:,end+1),f]=parametricEQ(filt.gain, filt.fc, filt.bandwidth,fs);
end

HA=H;

%%
load C_weight

EQs = uncalDev.Class1.C_weighted;

H=[];
for i=1:size(EQs.HPF)
    [H(:,end+1),f]=ResonantFilter(EQs.HPF(i).fc,'HPF',fs);
end
for i=1:size(EQs.LPF)
    [H(:,end+1),f]=ResonantFilter(EQs.LPF(i).fc,'LPF',fs);
end
for i=1:size(EQs.Parametric)
    filt=EQs.Parametric(i);
    [H(:,end+1),f]=parametricEQ(filt.gain, filt.fc, filt.bandwidth,fs);
end

HC = H;
%%

semilogx(A_weight(:,1), A_weight(:,2), C_weight(:,1), C_weight(:,2), f, sum(HA,2), f, sum(HC,2) )
axis([10 20000 -80 10])

legend('A-weight', 'C-weight', 'Approx. A-weight', 'Approx. C-weight');  

