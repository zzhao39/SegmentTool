I=imread('1.png');
I=rgb2gray(I);
BF=I;
strel_size1=15;

    se = strel('disk',strel_size1); %((((((((((((((param_to_be_tuned)))))))))))))))))))
    Ie = imerode(I, se);
    Iobr = imreconstruct(Ie, I);
    % figure, imshow(Iobr), title('Opening-by-reconstruction (Iobr)')
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    % figure, imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')

    % use Iobrcbr for marker extraction
    fgm = imregionalmax(Iobrcbr);
    % figure, imshow(fgm), title('Regional maxima of opening-closing by reconstruction(fgm)')
    I2 = I;
    I2(fgm) = 255;
    fgm4 = bwareaopen(fgm, strel_size1);
    I3 = I;
    I3(fgm4) = 255;
    % figure, imshow(I3), title('Modified regional maxima superimposed on original image (fgm4)')

    % Step 3: compute foreground markers
    bw = im2bw(Iobrcbr, graythresh(Iobrcbr));
    % figure, imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')

    D = bwdist(bw);
    DL = watershed(D);
    bgm = DL == 0;
    % figure, imshow(bgm), title('Watershed ridge lines (bgm)')

    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    % figure, imshow(gradmag,[]), title('Gradient magnitude (gradmag)')
    gradmag2 = imimposemin(gradmag, bgm | fgm4);
    % final watershed
    L = watershed(gradmag2);
    
    I_original=imread('1.png');
    
    background_label=L(1,1);
    I_r = I_original(:,:,1);
    I_g = I_original(:,:,2);
    I_b = I_original(:,:,3);
    [r,c] = size(I_b);


    for i=1:r
        for j=1:c
            if (I_b(i,j)>I_r(i,j)) && (I_b(i,j)>I_g(i,j)) 
                L(i,j) = background_label;
            end
        end
    end
    
Convex_threshold = 0.9;
    outImg=BF;
    
    folder_name = uigetdir;
    numRock = max(max(L));%the total number of segments 
    rockSize = [];%0*0 double array
    rockID = 1;
    %threshold = 0.91;
    negCount = 0;
    posCount = 0;
    [m,n] = size(outImg);%image BF is m*n matrix
    
    maxW = n/5; % upper bound
    maxH = maxW;
    L_color_label=L;
    for i = 1 : numRock
        wholeMask = L;
        [rows,cols] = find(L == i);
        wholeMask(find(L~=i)) = 0;
        A = max(min(rows)- 1, 1);
        B = min(max(rows)+ 1, m);
        C = max(min(cols)- 1, 1);
        D = min(max(cols)+ 1, n);

        height = max(rows) - min(rows)+ 1;
        width = max(cols)- min(cols) + 1;
        centerX = sum(rows)/size(rows,1);
        centerY = sum(cols)/size(cols,1);

        [K,area] = convhull(rows,cols);
        numPix = size(rows,1); %number of pixels;
        ratio = numPix/area;
        rockSize(rockID,1) = numPix;
        rockSize(rockID,2) = area; %the area of the convex hull
        rockSize(rockID,3) = numPix/area; %ratio
        singleMask = double(wholeMask(A:B, C:D));
        singleRock = double(outImg(A:B,C:D)).*double(singleMask);
        if ratio < Convex_threshold && height < maxH && width < maxW
            negCount = negCount + 1;
            rockSize(rockID,4) = 0; 
            imwrite(singleMask, [folder_name,'/FalseSeg',int2str(negCount),'_org_',int2str(rockID),'.jpg'],'jpg');
            L_color_label(find(L_color_label==i)) = 0;
        else
            posCount = posCount + 1;
            rockSize(rockID,4) = 1;
            imwrite(singleMask, [folder_name,'/TrueSeg',int2str(posCount),'_org_',int2str(rockID),'.jpg'],'jpg');   
        end
        rockID = rockID + 1;
    end
    Lrgb = label2rgb(L_color_label, 'jet', 'w', 'shuffle');
    figure, imshow(Lrgb), title('color labeled true/false segments')
    
    
    
    
    
a_lower=0.8;
a_upper=2.5;
b_lower=0.6;
b_upper=3;

ms = mrockSize(mrockSize(:,4)== 1,:);
r = ms(:,1)/ballSize/ballSize; % !!! make sure ballSize is updated for each input image

[x,y] = find(r >= a_lower & r <= a_upper);

p1 = sum(ms(x,1))/m/n;
[x1,y] = find(r < a_lower & r >= b_lower);
[x2,y] = find(r < b_upper & r >= a_upper);
p2 = (sum(ms(x1,1))+sum(ms(x2,1)))/m/n;