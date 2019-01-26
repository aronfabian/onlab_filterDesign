clear all
config()

Result_all = zeros(80,4);
for k = 1:4
    Result = zeros(80,1);
    switch k
            case 2
                fileNames = {'CutOf_Samsung Galaxy Ace 20170427_173930_SweepOnly_h_bSpectr'};
            case 3 
                fileNames = {'CutOf_Samsung Galaxy Alpha - 20170427_183001_SweepOnly_h_bSpectr'};
            case 4
                fileNames = {'CutOf_Samsung Galaxy Mini 2 - 20170427_174623_SweepOnly_h_bSpectr'};
    end
    for i = 1:4
        switch i
            case 2
                ERROR_SEL = 2;
            case 3 
                ERROR_SEL = 3;
            case 4
                ERROR_SEL = 4;
        end
        for j = 1:20
            [Filters,~] = invFilterDesign();
            Result((i-1)*20+j)=size(Filters,1);
            k
            i
            j
        end
    end
    Result_all(:,k)=Result;
end

save('filtnum.mat','Result_all')
%%
load('filtnum.mat','Result_all')
e1 = Result_all(1:20,:); % sum of abs (ERROR_SEL = 1)
e2 = Result_all(21:40,:); % sum of squares (ERROR_SEL = 2)
e3 = Result_all(41:60,:); % sum of abs where H_mic is outside the tolerance band (ERROR_SEL = 3)
e4 = Result_all(61:80,:); % sum of squares where H_mic is outside the tolerance band (ERROR_SEL = 4)

p1 = Result_all(:,1); % phone 1
p2 = Result_all(:,2); % phone 2
p3 = Result_all(:,3); % phone 3
p4 = Result_all(:,4); % phone 4

p1_e1 = p1(1:20);
p1_e2 = p1(21:40);
p1_e3 = p1(41:60);
p1_e4 = p1(61:80);
p2_e1 = p2(1:20);
p2_e2 = p2(21:40);
p2_e3 = p2(41:60);
p2_e4 = p2(61:80);
p3_e1 = p3(1:20);
p3_e2 = p3(21:40);
p3_e3 = p3(41:60);
p3_e4 = p3(61:80);
p4_e1 = p4(1:20);
p4_e2 = p4(21:40);
p4_e3 = p4(41:60);
p4_e4 = p4(61:80);

oszto1 = max(p1) - min(p1);
offset1 = min(p1) / oszto1;
oszto2 = max(p2) - min(p2);
offset2 = min(p2) / oszto2;
oszto3 = max(p3) - min(p3);
offset3 = min(p3) / oszto3;
oszto4 = max(p4) - min(p4);
offset4 = min(p4) / oszto4;

p1_norm = p1./oszto1 - offset1;
p2_norm = p2./oszto2 - offset2;
p3_norm = p3./oszto3 - offset3;
p4_norm = p4./oszto4 - offset4;

figure
hist(p1_norm(1:20))
hold 
hist(p1_norm(21:40))
hist(p1_norm(41:60))
hist(p1_norm(61:80))
figure
hist(p2_norm(1:20))
hold 
hist(p2_norm(21:40))
hist(p2_norm(41:60))
hist(p2_norm(61:80))
figure
hist(p3_norm(1:20))
hold 
hist(p3_norm(21:40))
hist(p3_norm(41:60))
hist(p3_norm(61:80))
figure
hist(p4_norm(1:20))
hold 
hist(p4_norm(21:40))
hist(p4_norm(41:60))
hist(p4_norm(61:80))
%%
 E1=[p1_norm(1:20); p2_norm(1:20); p3_norm(1:20); p4_norm(1:20)];
 E2=[p1_norm(21:40); p2_norm(21:40); p3_norm(21:40); p4_norm(21:40)];
 E3=[p1_norm(41:60); p2_norm(41:60); p3_norm(41:60); p4_norm(41:60)];
 E4=[p1_norm(61:80); p2_norm(61:80); p3_norm(61:80); p4_norm(61:80)];
 
 u1=unique(E1);
 u2=unique(E2);
 u3=unique(E3);
 u4=unique(E4);
 
 n1 = histc(E1,u1)/80;
 n2 = histc(E2,u2)/80;
 n3 = histc(E3,u3)/80;
 n4 = histc(E4,u4)/80;