% A-weight approximator

clear
load A_weight

fs=44100;
% g=7;
% [Ho,f]=parametricEQ(10,1000,1/3,44100);
paramEQs = [
     20, -44, 2.8
     28, -2.5,  0.8
     40, -8,  2.8
    100,  1,  2.8
   2000,  2,  2.8
   100000, -9.3, 0.45
    ];

H=[];
for i=1:size(paramEQs,1)
    [H(:,end+1),f]=parametricEQ(paramEQs(i,2),paramEQs(i,1),paramEQs(i,3),fs);
end


% fc=20; [H,f]=parametricEQ(-45,fc(end),2.8,fs);
% [H(:,end+1),f]=parametricEQ(-2,28,0.5,fs);
% [H(:,end+1),f]=parametricEQ(-8,40,2.8,fs);
% [H(:,end+1),f]=parametricEQ(1,100,2.8,fs);
% [H(:,end+1),f]=parametricEQ(2,2000,2.8,fs);
% [H(:,end+1),f]=parametricEQ(-9.3,100000,0.45,fs);
% [H3,f]=parametricEQ(g,1000*2^(1/9),1/9,44100);
% [H2,f]=parametricEQ(g,1000,1/6,44100);
semilogx(A_weight(:,1), A_weight(:,2), f, sum(H,2) )
axis([10 20000 -80 10])

