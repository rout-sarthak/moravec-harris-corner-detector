function [ moravec , harris ] = feature_detect( I )
%Function takes in input grayscale image and gives two pseudo colored
%using moravec and harris operator. The code was run on the three images
%given in previous homework assignments namely 'flower.pgm, tools.pgm,
%swan.pgm. It was also tested on a chequered pattern image. In all the
%cases performance oh harris operator was noticably better
figure(1); imshow(I);
%I = rgb2gray(I); %necessary to use if input image is not grayscale. But in
%the problem it is given that the image is grayscale
P = I;
Q = I;
numOfRows = size(I, 1);
numOfColumns = size(I, 2);
cornerness = zeros(numOfRows, numOfColumns);
[m,n]=size(I);
for i= 3 : (m-2)
    for j = 3: (n-2)
        threshold = 188; %threshold for cornerness
        A = I((i-1:i+1),(j-1:j+1));
        %shifting right
        R = I((i:i+2),(j-1:j+1));
        %difference of right square
        right_shift = abs(R-A).^2;
        %sum of right shift 
        SR = sum(sum(right_shift)); %sum function is used twice as the first sum gives
                                    %an array and whose's sum gives final answer

        %shifting left
        L = I((i-2:i),(j-1:j+1));
        %difference of left square
        left_shift = abs(L-A).^2;
        %sum of left shift 
        SL = sum(sum(left_shift));

        %shifting up
        U = I((i-1:i+1),(j-2:j));
        %difference of up square
        up_shift = abs(U-A).^2;
        %sum of up shift 
        SU = sum(sum(up_shift));

        %shifting down
        D = I((i-1:i+1),(j:j+2));
        %difference of down square
        down_shift = abs(D-A).^2;
        %sum of down shift 
        SD = sum(sum(down_shift));
        
        %shifting upleft
        UL = I((i-2:i),(j-2:j));
        %difference of up-left square
        upleft_shift = abs(UL-A).^2;
        %sum of upleft shift
        SUL = sum(sum(upleft_shift));
        
        %shifting upright
        UR = I((i:i+2),(j-2:j));
        %difference of up-right square
        upright_shift = abs(UR-A).^2;
        %sum of upright shift
        SUR = sum(sum(upright_shift));
        
        %shifting downleft
        DL = I((i-2:i),(j:j+2));
        %difference of down-left square
        downleft_shift = abs(DL-A).^2;
        %sum of downleft shift
        SDL = sum(sum(downleft_shift));

        %shifting downright
        DR = I((i:i+2),(j-2:j));
        %difference of down-right square
        downright_shift = abs(DR-A).^2;
        %sum of downright shift
        SDR = sum(sum(downright_shift));
        
        C = [SL, SR, SU, SD, SUL, SUR, SDL, SDR];
        cornerness(i,j) = min(C);
        if cornerness(i,j) < threshold
            cornerness(i,j) = 0;
        end   
    end
end

%Non-maximal Suppression
[Gmag,Gphase] = imgradient(cornerness);
[a,b]=size(cornerness);

for x = 2:a-1
    for y = 2:b-1
        O = Gphase(x,y);
        if O >= 157.5
            O = O - 180;
        end
        if O < -22.5
            O = O+180;
        end
        if (O >= -22.5) && (x < 22.5)
            neigh1 = Gmag(x-1,y);
            neigh2 = Gmag(x+1,y);
        else
            if (O >= 22.5) && (x < 67.5)
                neigh1 = Gmag(x-1,y-1);
                neigh2 = Gmag(x+1,y+1);
            else
                if (O >= 67.5) && (x < 112.5)
                    neigh1 = Gmag(x,y-1);
                    neigh2 = Gmag(x,y+1);
                else
                    if (O >= 112.5) && (x < 157.5)
                    neigh1 = Gmag(x-1,y+1);
                    neigh2 = Gmag(x+1,y-1);
                    end
                end
            end
        end
    if cornerness(x,y) >= neigh1 && cornerness(x,y)>= neigh2
        Glocalmax(x,y) = Gmag(x,y); 
    else
        Glocalmax(x,y) = 0;
    end
    if Glocalmax(x,y)>0 %all non zero points are corners
    P(x-2:x+2,y-2:y+2) = 255;
    Q(x-2:x+2,y-2:y+2) = 0;
    end
    end
end

redImage = Q;
greenImage = P;
blueImage = Q;
% Change the colors for the pixels that are lit
% Combine into a new RGB image.
moravec = cat(3, redImage, greenImage, blueImage);
figure(2); imshow(moravec);

%harris
[xx, yy] = meshgrid(-1:1, -1:1); %3x3 Window
sigma = 1;
Gxy = exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));

Gx = xx .* exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));
Gy = yy .* exp(-(xx .^ 2 + yy .^ 2) / (2 * sigma ^ 2));

P=I;
Q=I;
k = 0.04;

% image's x and y derivatives are taken
Ix = conv2(Gx, I);
Iy = conv2(Gy, I);
size(Ix);

% products of derivatives at every pixel is calculated here
Ix2 = Ix .^ 2;
Iy2 = Iy .^ 2;
Ixy = Ix .* Iy;

% sums of the products of derivatives at each pixel
Sx2 = conv2(Gxy, Ix2);
Sy2 = conv2(Gxy, Iy2);
Sxy = conv2(Gxy, Ixy);

Threshold = 50000;
numOfRows = size(I, 1);
numOfColumns = size(I, 2);
im = zeros(numOfRows, numOfColumns);
for x=1:numOfRows,
   for y=1:numOfColumns,
       x,y
       % matrix C
       C = [Sx2(x, y) Sxy(x, y); Sxy(x, y) Sy2(x, y)];
       
       % Harris Corner response
       HCR = det(C) - k * (trace(C) ^ 2);
       
       % Threshold on value of HCR
       if (HCR > Threshold)
          im(x, y) = HCR; 
       end
   end
end

%Non-maximal Suppression
[Gmag,Gphase] = imgradient(im);
[a,b]=size(im);
for x = 3:a-2
    for y = 3:b-2
        O = Gphase(x,y);
        if O >= 157.5
            O = O - 180;
        end
        if O < -22.5
            O = O+180;
        end
        if (O >= -22.5) && (x < 22.5)
            neigh1 = Gmag(x-1,y);
            neigh2 = Gmag(x+1,y);
        else
            if (O >= 22.5) && (x < 67.5)
                neigh1 = Gmag(x-1,y-1);
                neigh2 = Gmag(x+1,y+1);
            else
                if (O >= 67.5) && (x < 112.5)
                    neigh1 = Gmag(x,y-1);
                    neigh2 = Gmag(x,y+1);
                else
                    if (O >= 112.5) && (x < 157.5)
                    neigh1 = Gmag(x-1,y+1);
                    neigh2 = Gmag(x+1,y-1);
                    end
                end
            end
        end
    if im(x,y) >= neigh1 && im(x,y)>= neigh2
        Glocalmax(x,y) = Gmag(x,y); 
    else
        Glocalmax(x,y) = 0;
    end
    if Glocalmax(x,y)>0 %all non zero points are corners
    P(x-2:x+2,y-2:y+2) = 255;
    Q(x-2:x+2,y-2:y+2) = 0;
    end
    end
end

% Change the colors for the pixels that are lit
redImage = P;
greenImage = Q;
blueImage = Q;

% Combine into a new RGB image.
harris = cat(3, redImage, greenImage, blueImage);
figure(3); imshow(harris);
end
