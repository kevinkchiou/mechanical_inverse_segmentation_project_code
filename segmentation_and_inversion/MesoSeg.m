%{
Name1 = 'TestWing1.tif';
img = zeros([size(imread(Name1)),size(imfinfo(Name1),2)]);
SumIm = zeros(size(imread(Name1)));
for i = 1 : size(imfinfo(Name1),2)
    image = im2bw(imread(Name1,i));
    img(:,:,i) = bwareaopen(image,70);
    SumIm = SumIm+img(:,:,i);
end
%}
%%
%Name2 = 'MIPWing.tif';
Seg = struct('WSTil',[]);
for t = 1 : 1
    image=adapthisteq(imread('tiffs/apical_constrict_1_mod.tif'));h=fspecial('gaussian', 50,10);
    %image=adapthisteq(imread('tiffs/apical_constrict_2_mod.tif'));h=fspecial('gaussian',50,9);
    %image=adapthisteq(imread('tiffs/apical_constrict_3.tif'));h=fspecial('gaussian', 30,5);
    %h2 = fspecial('gaussian',40,30);
    im = imfilter(image,h);
    im2 = edge(image,'log',2/1000);
    seD = strel('disk',30);
    se0 = strel('line',2,0);
    se90 = strel('line',2,90);
    im2 = imerode(bwareaopen(imfill(imdilate(im2,seD),'holes'),10000),seD);
    %im2 = imerode(bwareaopen(imfill(imdilate(im2,seD),'holes'),200),seD);
    im2=imdilate(imdilate(imdilate(imerode(imerode(imerode(im2,seD),seD),seD),seD),seD),seD);
    im3 = im;
    im3(~im2) = Inf;
    BW = imregionalmin(im3);
    %BWcheck=BW;BWcheck(L==0)=255;figure(2);imshow(BWcheck);figure(1);
    %im2 = imimposemin(image,BW);
    %image = adapthisteq(image); %added in for better watershed hopefully
    L = watershed(imimposemin(image,BW));
    image2 = image;
    image2(L == 0) = 255;
    Lmod = watershed(imimposemin(image2,BW)); %done again since watershed finds edges but is leaky.
    %imshow(im2)
    %
    s = regionprops(L);
    Areas = cat(1,s.Area);
    Big = find(Areas > 1000);
    OK = find(Areas< 1000);
    for i = 1 : length(Big)
        L(L == Big(i)) = 1;
    end
    %imagesc(L); 
    Seg(t).WSTil = imclearborder(L);
end
figure(1);close(1);figure(1);imshow(image2);impixelinfo;
figure(2);close(2);figure(2);imshow(L);
%%
%{
v=getverts(L,1);
numcell=max(max(L));
gray2 = [gray;1 0 0];
%}
%%
%{
for t = 1 : 1
    image = imread(Name2,t); image=image(:,:,1); %weird reading of tif file
    L = Seg(t).WSTil;
   % imagesc(L);
     image(L == 0) = 256;
     L(bwlabel(imfill(imdilate(im2bw(256-image),[se90,se0]),'holes'))==0) = 1;
     image = imread(Name2,t);image=image(:,:,1); %weird reading of tif file
     image(L == 0) = 256;
     %imshow(image)
     hold on 
     s = regionprops(L);
     Centroids = cat(1,s.Centroid);
     plot(Centroids(:,1),Centroids(:,2),'r+')
     zoom(2);
     colormap(gray2);
%     camva(5);
    %view([0,90]);
    hold off
    M(t) = getframe;
end
%}
