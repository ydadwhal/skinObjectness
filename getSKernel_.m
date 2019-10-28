function [smKer]=getSKernel_(sigma)
% Get smooth filter
if sigma==0
    smKer=1;
    return;
end

dim=length(sigma); % row, column, third dimension
sz=max(ceil(sigma*2),1);
sigma=2*sigma.^2;

if dim==1
    d1=-sz(1):sz(1);
    
    smKer=exp(-((d1.^2)/sigma));
    
elseif dim==2
    [d2,d1]=meshgrid(-sz(2):sz(2),-sz(1):sz(1));
    
    smKer=exp(-((d1.^2)/sigma(1)+(d2.^2)/sigma(2)));
    
elseif dim==3
    [d2,d1,d3]=meshgrid(-sz(2):sz(2),-sz(1):sz(1),-sz(3):sz(3));
    
    smKer=exp(-((d1.^2)/sigma(1)+(d2.^2)/sigma(2)+(d3.^2)/sigma(3)));
    
else
    error('Not implemented');
end

smKer=smKer/sum(smKer(:));





