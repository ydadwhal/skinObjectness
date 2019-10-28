%this scrip tests the appliation of gabor filter on skin data whose output is fed to  salmat
%skin based gabor saliency
%Author : Yogeshwar Singh Dadwhal
%Date : 18-07-2019 (updated)
clear all;
[common] = Yinitialize();
fontSize=10
path=common.skindata;
flag=1;
% load ([path 'skinimageformdata.mat']);
load([ 'E:\MATLAB3\skinimageformdata.mat'])
wl=3;
% count=0;
count=20;

a=Pos;

sigma=1;

%        scales=4 and orientations=6
filterOutput = filterGabor(a(:,:,1),4,6,wl,2,0.65,1.5,1);
R=reshape([  real(filterOutput{1,1}),real(filterOutput{1,2}),real(filterOutput{1,3}),real(filterOutput{1,4}),real(filterOutput{1,5}),real(filterOutput{1,6}),real(filterOutput{2,1}),real(filterOutput{2,2}),real(filterOutput{2,3}),real(filterOutput{2,4}),real(filterOutput{2,5}),real(filterOutput{2,6}),real(filterOutput{3,1}),real(filterOutput{3,2}),real(filterOutput{3,3}),real(filterOutput{3,4}),real(filterOutput{3,5}),real(filterOutput{3,6}),real(filterOutput{4,1}),real(filterOutput{4,2}),real(filterOutput{4,3}),real(filterOutput{4,4}),real(filterOutput{4,5}),real(filterOutput{4,6}) ],1,24*size(a,1)*size(a,2));
filterOutput = filterGabor(a(:,:,2),4,6,wl,2,0.65,1.5,1);
G=reshape([  real(filterOutput{1,1}),real(filterOutput{1,2}),real(filterOutput{1,3}),real(filterOutput{1,4}),real(filterOutput{1,5}),real(filterOutput{1,6}),real(filterOutput{2,1}),real(filterOutput{2,2}),real(filterOutput{2,3}),real(filterOutput{2,4}),real(filterOutput{2,5}),real(filterOutput{2,6}),real(filterOutput{3,1}),real(filterOutput{3,2}),real(filterOutput{3,3}),real(filterOutput{3,4}),real(filterOutput{3,5}),real(filterOutput{3,6}),real(filterOutput{4,1}),real(filterOutput{4,2}),real(filterOutput{4,3}),real(filterOutput{4,4}),real(filterOutput{4,5}),real(filterOutput{4,6}) ],1,24*size(a,1)*size(a,2));
filterOutput = filterGabor(a(:,:,3),4,6,wl,2,0.65,1.5,1);
B=reshape([  real(filterOutput{1,1}),real(filterOutput{1,2}),real(filterOutput{1,3}),real(filterOutput{1,4}),real(filterOutput{1,5}),real(filterOutput{1,6}),real(filterOutput{2,1}),real(filterOutput{2,2}),real(filterOutput{2,3}),real(filterOutput{2,4}),real(filterOutput{2,5}),real(filterOutput{2,6}),real(filterOutput{3,1}),real(filterOutput{3,2}),real(filterOutput{3,3}),real(filterOutput{3,4}),real(filterOutput{3,5}),real(filterOutput{3,6}),real(filterOutput{4,1}),real(filterOutput{4,2}),real(filterOutput{4,3}),real(filterOutput{4,4}),real(filterOutput{4,5}),real(filterOutput{4,6}) ],1,24*size(a,1)*size(a,2));

M=cov([R;G;B]')
M
size(M);
[U,S,V] = svd(M);
Mtrans=U*(S^-1)*U'
clear ('filterOutput','R','G','B');
path=common.data;
% pat='C:\Users\YSDadwhal\Desktop\NEw Data\';
pat='F:\SDV2\4';
folderList1=dir([pat]);
folderList2=[folderList1.isdir];
folderList=folderList1(folderList2);
oro=size(folderList,1);
directory=1;%use this variable for batch processing for video
% for directory=1:15
    
    % path=([pat folderList(directory,1).name '\']);
    path=([pat '\']);
    
    filetype='*.jpg';
    fileList=dir([path filetype]);
    num=size(fileList,1);
    config.model=[];
    num=20;
    for i=num:num %change num to 1 for testing
        % ImgName=fileList(i,1).name
        ImgName=([num2str(i) '.jpg']);
        I=([path ImgName]);
        input_im=(double(imread((I))));
        input_im=imresize(input_im,0.4);
        input_im1=input_im;
        
        R=reshape(input_im1(:,:,1),1,size(input_im1,1)*size(input_im1,2));
        G=reshape(input_im1(:,:,2),1,size(input_im1,1)*size(input_im1,2));
        B=reshape(input_im1(:,:,3),1,size(input_im1,1)*size(input_im1,2));

        final=Mtrans*[R;G;B];
        final=reshape(final',size(input_im,1),size(input_im,2),size(input_im,3));
        [ skinSpatialMat ] = YogiSalLBP(final(:,:,1));
                  
       
        figure,
        subplot(2,2,1)
        imshow(final);
        fontSize=10
        title('Skin Distinctiveness', 'FontSize', fontSize);
        % level = graythresh(skinSpatialMat)
%         binaryImage=im2bw(skinSpatialMat,0.4);
          binaryImage=imbinarize(skinSpatialMat,'global');
        subplot(2,2,2);
        % imshow(binaryImage);
        % title('Binary Image', 'FontSize', fontSize);
        imshow(skinSpatialMat);
        title('Spatial Saliency', 'FontSize', fontSize);
        % B = bwboundaries(binaryImage);
        % imshow(binaryImage)
        
        binaryImage = imfill(binaryImage, 'holes');
        
        % Remove tiny regions.
        
%       [binaryImage config]=objectnessTracking(binaryImage,config);
        
        
        binaryImage = bwareaopen(binaryImage, 100);% set mimimun number of pizels to be deleted as blob
      
%       imwrite((binaryImage),['H:\DATASET\VICTIM data\MOHALI VICTIM VIDEO DATA WITH GT\' num2str(directory) '\BinTrackGab1\' ImgName ]);
        
        % Extract the largest area using ImageAnalyst's custom function ExtractNLargestBlobs().
        % biggestBlob=binaryImage;
        biggestBlob = ExtractNLargestBlobs(binaryImage, 4);%number of blobs to be extracted
        % Display the image.
        % subplot(2,2,4);
        % imshow(biggestBlob);
        % title('Final Image', 'FontSize', fontSize);
        
        subplot(2,2,3);
        imshow(binaryImage);
        str = sprintf('Binary Image Folder=%d Number= %s',directory,ImgName);
        title(str,'FontSize',fontSize);
        %
        % subplot(2,3,5);
        % imshow(biggestBlob);
        % title('Post processed Binary Image','FontSize',fontSize);
        
        
        subplot(2,2,4);
        imshow(uint8(input_im));
        title('Skin Objectness','FontSize',fontSize);
        
        
        [labeledImage, numberOfBlobs] = bwlabel(biggestBlob, 8);
        % Get all the blob properties.
        blobMeasurements = regionprops(labeledImage, 'BoundingBox','Area');
        allBlobAreas = [blobMeasurements.Area];
        % Display the original gray scale image.
        subplot(2,2,4);
        % Loop through all blobs, putting up Bounding Box.
        hold on; % Prevent boxes from blowing away the image and prior boxes.
        
        for k = 1 : numberOfBlobs
            k
            boundingBox = blobMeasurements(k).BoundingBox;   % Get box.
            x1 = boundingBox(1)-4;
            y1 = boundingBox(2)-4;
            x2 = x1 + boundingBox(3) - 1+4;
            y2 = y1 + boundingBox(4) - 1+4;
            verticesX = [x1 x2 x2 x1 x1];
            verticesY = [y1 y1 y2 y2 y1];
            plot(verticesX, verticesY, 'r-', 'LineWidth', 1);
        end
        hold off;
        
        
        % str='C:\Users\Yogeshwar\Desktop\saliency data\Cheng and Mitra Saliency\GT\MSRA10K_Imgs_GT\GT\';
        % im=ImgName((strfind(ImgName,'_')+1):end-3);
        % img=([str im 'png']);
        % gt=imread(img);
        % subplot(2,3,3);
        % imshow(uint8(gt));
        % title('Ground Truth','FontSize',fontSize);
        
        
        % subplot(2,3,6);
        % imshow(uint8(input_im));
        % title('GT bounding box','FontSize',fontSize);
        
        % [labeledImage, numberOfBlobs] = bwlabel(gt, 8);
        % % Get all the blob properties.
        % blobMeasurements = regionprops(labeledImage, 'BoundingBox','Area');
        % allBlobAreas = [blobMeasurements.Area];
        % % Display the original gray scale image.
        % subplot(2,3,6);
        % % Loop through all blobs, putting up Bounding Box.
        % hold on; % Prevent boxes from blowing away the image and prior boxes.
        %
        % for k = 1 : numberOfBlobs
        % k
        %     boundingBox = blobMeasurements(k).BoundingBox;   % Get box.
        %       x1 = boundingBox(1);
        %       y1 = boundingBox(2);
        %       x2 = x1 + boundingBox(3) - 1;
        %       y2 = y1 + boundingBox(4) - 1;
        %       verticesX = [x1 x2 x2 x1 x1];
        %       verticesY = [y1 y1 y2 y2 y1];
        %       plot(verticesX, verticesY, 'r-', 'LineWidth', 2);
        % end
        % hold off;
        
        
        
        set(gcf, 'visible','off')
        % ImgName=[ImgName(1:end-3) 'jpg'];
        pa=[path 'test1\' ImgName ] %skindistict Sk
        saveas(gcf,pa);
        close all;
        imshow(uint8(input_im));
        hold on;
        for k = 1 : numberOfBlobs
            k
            boundingBox = blobMeasurements(k).BoundingBox;   % Get box.
            x1 = boundingBox(1);
            y1 = boundingBox(2);
            x2 = x1 + boundingBox(3) - 1;
            y2 = y1 + boundingBox(4) - 1;
            verticesX = [x1 x2 x2 x1 x1];
            verticesY = [y1 y1 y2 y2 y1];
            plot(verticesX, verticesY, 'r-', 'LineWidth', 1);
        end
        hold off
        set(gcf, 'visible','off')
        pa=[path 'test1\' 'SD\' ImgName ]; %skindistict Sk
        saveas(gcf,pa);
        
        clear ImgName R S U V W a binaryImage binaryImageB binaryImageR binaryImageG B filterOutput G I M final fontSize level maxit Pos Neg;
        clear common allBlobAreas biggestBlob blobMeasurements boundingBox final input_im
        clear labeledImage maxit numberOfBlobs skinSpatialMat
        % uiwait(msgbox('Done with Yogeshwar'));
        close all;
        
    end
    
    
    
% end

% end
