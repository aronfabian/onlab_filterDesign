param(1).gain = -20.62;
param(1).fc = 50.0;
param(1).bandwidth = 3.9656; 

param(2).gain = -17.90;
param(2).fc = 15600.0;
param(2).bandwidth = 1.2894; 

% Ezzel lehet majd 1 utasításban az param struktúra tömböt létrehozni: 
% foo = struct('field_a',num2cell([1,2,3,4]), 'field_b',num2cell([4,8,12,16]))

hpf(1).fc = 100;
hpf(2).fc = 50;

lpf(1).fc = 15000;

filters.HPF = hpf;
filters.LPF = lpf;
filters.Parametric = param;

Class1.A_weighted = filters;
Class1.C_weighted = filters;
Class2.A_weighted = filters;
Class2.C_weighted = filters;

Phone(1).model = '12345';
Phone(1).market_name = 'OnePlus3T';
Phone(1).phone_id = '67890'; % Android ID
Phone(1).imei = '123456789012345'; % 
Phone(1).Class1 = Class1 
Phone(1).Class2 = Class2 

Phone(2).model = 'SM-A320FL';
Phone(2).market_name = 'Samsung Galaxy A3(2017)';
Phone(2).phone_id = '67890';
Phone(2).imei = '123456789012345'; % 
Phone(2).Class1 = Class1 
Phone(2).Class2 = Class2 

% Phone1.androidID = '12345'
% Phone1.type = 'OnePlus3T';
% Phone1.A_weighted = filters;
% 
% Phone2.androidID = '67890'
% Phone2.type = 'Galaxy S8';
% Phone2.C_weighted = filters;

% PhoneDatabase.database = Phone

% json_txt=savejson('Phone',Phone); % Ez nem menti a fájlt, csak visszatér a txt-vel
json_txt=savejson('Phones',Phone,'PhoneDatabase.json'); % Ez nem menti a fájlt, csak visszatér a txt-vel

