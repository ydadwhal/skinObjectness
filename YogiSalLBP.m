function [ skinSpatialMat ] = YogiSalLBP( I )
%Author: Yogeshwar Singh Dadwhal
%Date: 19-04-2016

I=padarray(I,[25,25],'replicate');
if size(I,3)>1
Igray=rgb2gray(I);
else
Igray=I;
end
[aaa bbb]=size(Igray);
% Igray=zeros(size(Igray));
 
% Igray = double(YogiLBP(Igray));

%relative window size for multiple scales
config.windowRows=[0.25 0.3 0.5 0.7]; 
config.windowCols=[0.1 0.3 0.5 0.4]; 
config.sampleStep=[0.01 0.015 0.03 0.04];
    
config.verticalBorder=0.1; 
config.horizontalBorder=0.1;
config.pH1=0.2; 

feat(1).numBin=60; 
feat(1).minVal=0;  % Minimum feature value
feat(1).maxVal=ceil(max(max(Igray))); % Maximum feature value
feat(1).smoothFactor=5; % Smoothing factor
% Other
config.verbose=2; 
config.flipBox=0; 
config.fullImage=[0 0 0 2]; % 0-> whole window inside image, 1-> kernel inside image, 2->full image. Corresponds to each feature in dataCell 


% Initialize
config.numRow=size(Igray,1);
config.numCol=size(Igray,2);
config.numFeat=3;
if config.numCol>config.numRow && config.flipBox==1
    tmp=config.windowRows; config.windowRows=config.windowCols; config.windowCols=tmp;
end   

numWindow=length(config.windowRows);

% Window rows and columns
config.windowRows=round(config.windowRows*max(config.numRow,config.numCol));
config.windowCols=round(config.windowCols*max(config.numRow,config.numCol));

% Window step size
if(length(config.sampleStep)<numWindow)
    config.sampleStep=repmat(config.sampleStep,1,numWindow);
end
config.sampleStep=ceil(config.sampleStep*max(config.numRow,config.numCol));

% Border lengths
if(length(config.verticalBorder)<numWindow)
    config.verticalBorder=repmat(config.verticalBorder,1,numWindow);
end
if(length(config.horizontalBorder)<numWindow)
    config.horizontalBorder=repmat(config.horizontalBorder,1,numWindow);
end
config.verticalBorder=round(config.verticalBorder.*max(config.windowRows,config.windowCols));
config.horizontalBorder=round(config.horizontalBorder.*max(config.windowRows,config.windowCols));

% Increase window sizes by border sizes
config.windowRows=config.windowRows+2*config.verticalBorder;
config.windowCols=config.windowCols+2*config.horizontalBorder;

% Full image flag
if(length(config.fullImage)<numWindow)
    config.fullImage=repmat(config.fullImage(1),[1,numWindow]);
end


numFeat=3;
numRow=size(Igray,1);
numCol=size(Igray,2);
skinSpatialMat=zeros(numRow,numCol,length(config.windowRows));%fonfig.windows rows are 4
intHist=cell(1,numFeat);
smoothingKernel=cell(1,numFeat);
finalMat=cell(1,numFeat);
if nargout<2
    config.compIndSal=0;
else
    config.compIndSal=1;
end
% Compute integral histogram image
    [intHist{1},Igray]=formIntegralHist_(Igray,feat(1).numBin,feat(1).minVal,feat(1).maxVal);
    smoothingKernel{1}=getSKernel_(feat(1).smoothFactor);
    if config.compIndSal>0
        finalMat{1}=zeros(numRow,numCol);
    end

   inthistsize=size(intHist)
   
%  sliding windows through image
for winSizeId=1:length(config.windowRows) % 1:4
    % Initialize window parameters
    halfRow=floor(config.windowRows(winSizeId)/2);
    halfCol=floor(config.windowCols(winSizeId)/2);
    halfKerRow=halfRow-config.verticalBorder(winSizeId);
    halfKerCol=halfCol-config.horizontalBorder(winSizeId);
    kernelRowIdx=-halfKerRow:halfKerRow;
    kernelColIdx=-halfKerCol:halfKerCol;
    
    % Initialize window positions
    if sum(config.fullImage(1:numFeat))>0
        rowCent=(config.sampleStep(winSizeId)+1):config.sampleStep(winSizeId):numRow;
        colCent=(config.sampleStep(winSizeId)+1):config.sampleStep(winSizeId):numCol;
        rowCent(rowCent>(numRow))=[];
        colCent(colCent>(numCol))=[];
        numWindow=length(rowCent)*length(colCent);
    else
        rowCent=(halfRow+1):config.sampleStep(winSizeId):numRow;
        colCent=(halfCol+1):config.sampleStep(winSizeId):numCol;
        rowCent(rowCent>(numRow-halfRow))=[];
        colCent(colCent>(numCol-halfCol))=[];
        numWindow=length(rowCent)*length(colCent);
    end

    pH1=config.pH1;%0.2
    pH0=1-pH1;%0.8

    cnt=1;
    for winRowPos=rowCent
        for winColPos=colCent
            
            cnt=cnt+1;
                  
            k11=[winRowPos-halfKerRow-1,winColPos-halfKerCol-1]+1; 
            k12=[winRowPos-halfKerRow-1,winColPos+halfKerCol]+1;
            k21=[winRowPos+halfKerRow,winColPos-halfKerCol-1]+1;
            k22=[winRowPos+halfKerRow,winColPos+halfKerCol]+1;
            w11=[winRowPos-halfRow-1,winColPos-halfCol-1]+1;
            w12=[winRowPos-halfRow-1,winColPos+halfCol]+1;
            w21=[winRowPos+halfRow,winColPos-halfCol-1]+1;
            w22=[winRowPos+halfRow,winColPos+halfCol]+1;
            
            if sum([k11<1 k12(1)<1 k12(2)>(numCol+1) k21(1)>(numRow+1) k21(2)<1 k22(1)>(numRow+1) k22(2)>(numCol+1)])>0
                innerWindowFlag=2;
            elseif sum([w11<1 w12(1)<1 w12(2)>(numCol+1) w21(1)>(numRow+1) w21(2)<1 w22(1)>(numRow+1) w22(2)>(numCol+1)])>0
                innerWindowFlag=1;
            else
                innerWindowFlag=0;
            end
            
            k11=min(max(k11,[1,1]),[numRow+1,numCol+1]);
            k12=min(max(k12,[1,1]),[numRow+1,numCol+1]);
            k21=min(max(k21,[1,1]),[numRow+1,numCol+1]);
            k22=min(max(k22,[1,1]),[numRow+1,numCol+1]);
            
            w11=min(max(w11,[1,1]),[numRow+1,numCol+1]);
            w12=min(max(w12,[1,1]),[numRow+1,numCol+1]);
            w21=min(max(w21,[1,1]),[numRow+1,numCol+1]);
            w22=min(max(w22,[1,1]),[numRow+1,numCol+1]);
            
            featureRowIndex=kernelRowIdx+winRowPos;
            featureRowIndex(featureRowIndex<1 | featureRowIndex>numRow)=[];
            featureColIndex=kernelColIdx+winColPos;
            featureColIndex(featureColIndex<1 | featureColIndex>numCol)=[];
            
            % Calculate histograms
            PrH1=pH1;
            PrH0=pH0;
            
            %yogi for Igray starts
            for i=1:1
                if innerWindowFlag==0 || config.fullImage(i)>=innerWindowFlag
                    pX_H1=double(intHist{i}(k22(1),k22(2),:)+intHist{i}(k11(1),k11(2),:)-intHist{i}(k12(1),k12(2),:)-intHist{i}(k21(1),k21(2),:));
                    pX_H0=double(intHist{i}(w22(1),w22(2),:)+intHist{i}(w11(1),w11(2),:)-intHist{i}(w12(1),w12(2),:)-intHist{i}(w21(1),w21(2),:))-pX_H1;
               pX_H1Size= size(pX_H1)
                    pX_H1=filterDistribution_(smoothingKernel{i},pX_H1(:)',feat(i).numBin);
                    pX_H0=filterDistribution_(smoothingKernel{i},pX_H0(:)',feat(i).numBin);
                
                    PrX_H1=pX_H1(Igray(featureRowIndex,featureColIndex));
                    PrX_H0=pX_H0(Igray(featureRowIndex,featureColIndex));
                    
                    PrX_H1Size=size(PrX_H1)
                    
                    PrH1=PrH1.*PrX_H1;
                    PrH0=PrH0.*PrX_H0;
                    
                    PrH1Size=size(PrH1)
                
                    if config.compIndSal>0
                        PrH1_X_perFeat=(PrX_H1*pH1)./(PrX_H1*pH1+PrX_H0*pH0);
                        finalMat{i}(featureRowIndex,featureColIndex)=max(PrH1_X_perFeat,finalMat{i}(featureRowIndex,featureColIndex));
                    end
                end
            end
             
            
            
            PrH1_X=PrH1./(PrH1+PrH0);
            
            skinSpatialMat(featureRowIndex,featureColIndex,winSizeId)=max(PrH1_X,skinSpatialMat(featureRowIndex,featureColIndex,winSizeId));

        end
    end
    if(config.verbose>=1)
        fprintf('\n');
    end
end

%%%%%YOGI for 


size(skinSpatialMat)
% Concatenate saliency images w;ith different resolutions
skinSpatialMat=max(skinSpatialMat,[],3);
skinSpatialMat=skinSpatialMat(1+25:end-25,1+25:end-25);
% 

% p=waitforbuttonpress()
end
% end



