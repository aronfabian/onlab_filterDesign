%% Generate and plot responses
fs = 44100; 

% Speaker response 
[Ho,f]=parametricEQ(6,1000,2,fs);
[H1,f]=parametricEQ(-4,7000,1/2,fs);
[H3,f]=parametricEQ(4.5,8000,1/2,fs);
[H2,f]=parametricEQ(-5,16000,1/2,fs);
[H4,f]=parametricEQ(1,200,1,fs);
[H5,f]=parametricEQ(2,20000,1,fs);

H_spk = Ho+H1+H2+H3+H4+H5;

% Target response
[H_trgt,f]=parametricEQ(0,1000,100,fs);

% plot
figure(1)
semilogx(f, [H_spk H_trgt])
axis([10 30000 -10 10])
title('Speaker and target responses')
legend('H-speaker', 'H-target')
xlabel('frequency [Hz]')
ylabel('gain [dB]')

%% Calc

% Error areas
cross = zeros(size(f,1),1);
areaNum = 1;
for i = 2:size(f,1)
    if ((H_spk(i) > H_trgt(i)) && (H_spk(i-1) < H_trgt(i-1)))
        cross(i)=1;
        areaNum = areaNum + 1;
    end
    if ((H_spk(i) < H_trgt(i)) && (H_spk(i-1) > H_trgt(i-1)))
        cross(i)=1;
        areaNum = areaNum + 1;
    end
end

errorAreas = zeros(areaNum,1);

m = 1;
for i = 1:size(f,1)
    if (cross(i) == 1)
        m = m + 1;
    end
    errorAreas(m) = errorAreas(m) + (H_spk(i) - H_trgt(i));
end