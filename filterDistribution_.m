function dist=filterDistribution_(filterKernel,dist,numBin)
% Smooth distribution
if numel(filterKernel)==1
    dist=dist(:)/sum(dist(:));
    return;
end

numDim=length(numBin);

if numDim==1
    %smoothDist=conv(dist,filterKernel,'same');
    
    lenDist=length(dist);
    hlenKernel=(length(filterKernel)-1)/2;

    dist=[dist(1)*ones(1,hlenKernel),dist,dist(end)*ones(1,hlenKernel)];
    dist=conv(dist,filterKernel);
    lenSmoothDist=length(dist);
    offset=(lenSmoothDist-lenDist)/2;
    dist=dist((offset+1):(lenSmoothDist-offset));
    
elseif numDim==2
    dist=reshape(dist,numBin);
    
    dist=conv2(filterKernel,filterKernel,dist,'same');
    
else
    dist=reshape(dist,numBin);
    
    for i=1:numDim
        fker=ones(1,numDim);
        fker(i)=length(filterKernel);
        fker=zeros(fker);
        fker(:)=filterKernel(:);
        
        dist=convn(dist,fker,'same');
        
    end
    
end

dist=dist(:)/sum(dist(:));






