%% agglomerative hierarchical cluster analysis AHCA
%% Refenence
% M. Miljkovi?, T. Chernenko, M. J. Romeo, B. Bird, C. Matthäus, and M. Diem, "Label-free imaging of human cells: algorithms for image 
% reconstruction of Raman hyperspectral datasets," Analyst 135, 2002 (2010).

%  https://www.mathworks.com/help/stats/clusterdata.html
%  https://www.mathworks.com/help/stats/pdist.html
%  https://www.mathworks.com/help/stats/hierarchical-clustering.html
clear,clc,close all
clst_n = 5; %bkg, cytoplasma, nucleus, ER, lipid
img_name ='001.tif';%read the hyperspectral stack 
img = read_tiff(img_name);
img = img(1:200,1:200,:);
size_x = size(img,1);
size_y = size(img,2);
lambda_n = size(img,3);
spct_n = size_x * size_y;
img_2d = reshape(img,[size_x*size_y,lambda_n]);
%%
T1 = clusterdata(img_2d,'Cutoff',1);
img_clsted = reshape(T1,[size_x,size_y]);
imshow(img_clsted)
%%
img_color = psudo_clr(img_clsted,clst_n);           
img_color_re = reshape(img_color,[size_x,size_y,3]);
imshow(img_color_re)
%% test pdist, linkage and dendrogram
% how to determine cutoff
Y = pdist(img_2d);
Z = linkage(Y);
dendrogram(Z);
T = cluster(Z,'maxclust',clst_n);%'cutoff',0.58);%
img_clsted = reshape(T,[size_x,size_y]);
figure,imshow(img_clsted)
%%
function img_stack = read_tiff(img_name)
    info = imfinfo(img_name);
    num_imgs = length(info);
    img_stack = [];
    for kk = 1:num_imgs
        currentimg = imread(img_name,kk,'Info',info);
        img_stack(:,:,kk) = currentimg;
    end
end
% function img_color = psudo_clr(img_clsted,clst_n)
%     color_code = [0,0,0;0,0,1;0,1,1;0,1,0;1,0,0];
%     img_color = zeros([size(img_clsted(:),1),3]);
%     for kk = 1:clst_n
%         locs = find(img_clsted(:) == kk);
%         img_color(locs,:) = repmat(color_code(kk,:),length(locs),1);%color_code(kk,:);
%     end
% end

