function L = getCost2( Image )
%%%%%%%%% Calculate local cost of 8 links %%%%%%%%

%LINK 1
ker(:,:,1)=1/4 .* [0 1 1; 0 0 0 ; 0 -1 -1];
%LINK 2
ker(:,:,2)=1/sqrt(2) .* [0 -1 0; 0 0 1; 0 0 0];
%LINK 3
ker(:,:,3)=1/4 .* [1 0 -1; 1 0 -1; 0 0 0];
%LINK 4
ker(:,:,4)=1/sqrt(2) .* [0 1 0; -1 0 0; 0 0 0];
%LINK 5
ker(:,:,5)=1/4 .* [1 1 0; 0 0 0; -1 -1 0];
%LINK 6
ker(:,:,6)=1/sqrt(2) .* [0 0 0; 1 0 0; 0 -1 0];
%LINK 7
ker(:,:,7)=1/4 .* [0 0 0; 1 0 -1; 1 0 -1];
%LINK 8
ker(:,:,8)=1/sqrt(2) .* [0 0 0; 0 0 -1; 0 1 0];

DLink = zeros(size(Image,1),size(Image,2),8);

if size(Image,3) == 3 %RGB image
    red = double(Image(:,:,1));
    green = double(Image(:,:,2));
    blue = double(Image(:,:,3));

    result = zeros(size(Image,1),size(Image,2),8,3);
    for i = 1:8
        result(:,:,i,1) = abs(imfilter(red,ker(:,:,i)));
        result(:,:,i,2) = abs(imfilter(green,ker(:,:,i)));
        result(:,:,i,3) = abs(imfilter(blue,ker(:,:,i)));
    end

    DLink = sqrt((result(:,:,:,1) .* result(:,:,:,1) +...
        result(:,:,:,2) .* result(:,:,:,2) +...
        result(:,:,:,3) .* result(:,:,:,3))./3);
    
else size(Image,3) == 1 %grayscale & others
    for i = 1:8
        DLink(:,:,i,1) = abs(imfilter(Image,ker(:,:,i)));
    end
end

maxD = max(max(max(DLink)));
L = maxD - DLink;
L(:,:,2) = sqrt(2)*L(:,:,2);
L(:,:,4) = sqrt(2)*L(:,:,4);
L(:,:,6) = sqrt(2)*L(:,:,6);
L(:,:,8) = sqrt(2)*L(:,:,8);
% L is the same dimension as DLink

end

