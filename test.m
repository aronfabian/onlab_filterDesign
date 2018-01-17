clear all
clc

 fileIDW = fopen('results.txt','w');fclose(fileIDW);
for i = 1:4
     text = fileread('results.txt');
     fileIDW = fopen('results.txt','w');
    switch(i)
        case 1
            hpf = 0;
            lpf = 0;
            fprintf(fileIDW,[text,'\n*********************************************************************************\n*********************************************************************************','\nHPF, LPF = 0\n']);
        case 2 
            hpf = 1;
            lpf = 0;
            fprintf(fileIDW,[text,'\n*********************************************************************************\n*********************************************************************************','\nHPF = 1, LPF = 0\n']);
        case 3
            hpf = 0;
            lpf = 1;
            fprintf(fileIDW,[text,'\n*********************************************************************************\n*********************************************************************************','\nHPF = 0, LPF = 1\n']);
        case 4
            hpf = 1;
            lpf = 1;

            fprintf(fileIDW,[text,'\n*********************************************************************************\n*********************************************************************************','\nHPF, LPF = 1\n']);
    end
    fclose(fileIDW);
    
    
    configtest(hpf,lpf)
    
    for n = 1:10
        invFilterDesign;
        text1 = fileread('output.txt');
        text2 = fileread('results.txt');
        fileIDW = fopen('results.txt','w');
        fprintf(fileIDW,[text2,'\n**\n',text1]);
        fclose(fileIDW);
        
    end
end