function output = filterGabor(im, scales, orientation, minWL, mult,sigmaf, dTheta, returnVal)
   
[rows cols] = size(im);					
imFFT = fft2(im);                 
output = cell(scales, orientation);         

[x,y] = meshgrid( [-cols/2:(cols/2-1)]/cols,...
		  [-rows/2:(rows/2-1)]/rows);
dis = sqrt(x.^2 + y.^2);                  
dis(round(rows/2+1),round(cols/2+1)) = 1; % Get rid of the 0 dis value in the middle 

theta = atan2(-y,x);             
sintheta = sin(theta);
costheta = cos(theta);
clear x; clear y; clear theta;      
thetaSigma = pi/orientation/dTheta;  

for o = 1:orientation,                  
  angl = (o-1)*pi/orientation;      %Filter angle     
  wavelength = minWL;        % filter wavelength.

  ds = sintheta * cos(angl) - costheta * sin(angl);     % Difference in sine.
  dc = costheta * cos(angl) + sintheta * sin(angl);     % Difference in cosine.
  dtheta = abs(atan2(ds,dc));                           % Absolute angular distance.
  spread = exp((-dtheta.^2) / (2 * thetaSigma^2));      

  for s = 1:scales,                  
    fo = 1.0/wavelength;                  % Centre frequency of filter.
    logGabor = exp((-(log(dis/fo)).^2) / (2 * log(sigmaf)^2));  
    logGabor(round(rows/2+1),round(cols/2+1)) = 0; 
% size(logGabor)
% sp=size(spread)
    filter = fftshift(logGabor .* spread);
    output{s,o} = ifft2(imFFT .* filter);    
    wavelength = wavelength * mult;       
  end                                    

end  



