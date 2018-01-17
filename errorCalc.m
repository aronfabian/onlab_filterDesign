function [errorArea,errorTol] = errorCalc(H_mic,H_trgt,tol_interp,configFile)
%ERRORCALC Calculates difference between two vectors
%   It can calculate:
%                       sum of abs (sumConf = 1)
%                       sum of squares (sumConf = 2)
%                       sum of abs where H_mic is outside the tolerance band (sumConf = 3)
%                       sum of squares where H_mic is outside the tolerance band (sumConf = 4)

load(configFile)

switch(ERROR_SEL)
    case 1
        errorArea = sum(abs(H_mic-H_trgt));
    case 2
        errorArea = sum((H_mic-H_trgt).^2);
    case 3
        errorArea = sum(H_mic((H_mic > tol_interp(1,:)) | (H_mic < tol_interp(2,:))) - H_trgt((H_mic > tol_interp(1,:)) | (H_mic < tol_interp(2,:)))); 
    case 4
        errorArea = sum((H_mic((H_mic > tol_interp(1,:)) | (H_mic < tol_interp(2,:))) - H_trgt((H_mic > tol_interp(1,:)) | (H_mic < tol_interp(2,:)))).^2); 
end
errorTol = length(find((H_mic > tol_interp(1,:)) | (H_mic < tol_interp(2,:))));
end

