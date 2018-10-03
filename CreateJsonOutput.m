function [] = CreateJsonOutput(filters)

p_i = 0;
h_i = 0;
l_i = 0; 
for i = 1 : size(filters,1)
    switch(filters(i,6))
        case 0 % parametric
            p_i = p_i + 1;
            param(p_i).gain = filters(i,1);
            param(p_i).fc = filters(i,2);
            param(p_i).bandwidth = filters(i,3);
        case 1 % hpf
            h_i = h_i + 1;
            hpf(h_i).fc = filters(i,2);
            
        case 2 % lpf
            l_i = l_i + 1;
            lpf(l_i).fc = filters(i,2);
    end
    
end
filt.HPF = hpf;
filt.LPF = lpf;
filt.Parametric = param;
Phone(1).androidID = '12412411'
Phone(1).type = 'OnePlus3T';
Phone(1).A_weighted = filt;
% param(1).gain = -20.62;
% param(1).fc = 50.0;
% param(1).bandwidth = 3.9656; 
% 
% param(2).gain = -17.90;
% param(2).fc = 15600.0;
% param(2).bandwidth = 1.2894;
% 
% % Ezzel lehet majd 1 utasításban az param struktúra tömböt létrehozni:
% % foo = struct('field_a',num2cell([1,2,3,4]), 'field_b',num2cell([4,8,12,16]))
% 
% hpf(1).fc = 100;
% hpf(2).fc = 50;
% 
% filters.HPF = hpf;
% filters.Parametric = param;
% 
% 
% Phone(1).androidID = '12345'
% Phone(1).type = 'OnePlus3T';
% Phone(1).A_weighted = filters;
% 
% Phone(2).androidID = '67890'
% Phone(2).type = 'Galaxy S8';
% Phone(2).A_weighted = filters;
% 
% 
% PhoneDatabase.database = Phone
% 
% json_txt=savejson('Phones',Phone,'PhoneDatabase.json'); 
% 
