function [ out config ] = yogiTracking( binaryImage,config )

% [labeledImage, numberOfBlobs] = bwlabel(binaryImage, 8);
% blobMeasurements = regionprops(labeledImage, 'BoundingBox','Area','FilledArea');
% boundingBox = blobMeasurements.BoundingBox;
% 
%     if isempty(model)
%         model=binaryImage;
%     else
%         model=not(xor((model),(binaryImage)));
%     end
% out=(model)&(binaryImage);

config.k=5; %memory of last k frames 
config.windowRows=[0.25 0.3 0.5 0.7]; % Window row size relative to image size max(imRow,imCol))
config.windowCols=[0.1 0.3 0.5 0.4]; % Window column size relative to image size max(imRow,imCol))
config.sampleStep=[0.01 0.015 0.03 0.04];
config.alpha=0.5;

    if isempty(config.model)
        [m n]=size(binaryImage);
        config.model=zeros(m,n,config.k);
%         config.model=zeros(m,n);
        config.ctr=0;
    end
config.ctr=config.ctr+1;
if config.ctr>5
    config.ctr=1;
end

%model(:,:,config.ctr)=binaryImage;
% model=config.model
% M=not(xor((model),(binaryImage)));
% model=config.alpha*M+(1-config.alpha)*model;
% out=(model)&(binaryImage);
% config.model=model;


model=config.model
% M=not(xor((model),(binaryImage)));
% model=config.alpha*M+(1-config.alpha)*model;
model(:,:,config.ctr)=binaryImage;
out=model(:,:,1)+model(:,:,2)+model(:,:,3)+model(:,:,4)+model(:,:,5);
out=uint8(round((out-2)./3));
config.model=model;


end

