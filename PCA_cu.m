%% PCA
%% Refenence
% M. Miljkovi?, T. Chernenko, M. J. Romeo, B. Bird, C. Matthäus, and M. Diem, "Label-free imaging of human cells: algorithms for image 
% reconstruction of Raman hyperspectral datasets," Analyst 135, 2002 (2010).

% clear,clc,close all
% clst_n = 5; %bkg, cytoplasma, nucleus, ER, lipid
% img_name ='001.tif';%read the hyperspectral stack 
% img = read_tiff(img_name);
% img = img(1:200,1:200,:);
function img_clsted = PCA_cu(img,clst_n)
size_x = size(img,1);
size_y = size(img,2);
lambda_n = size(img,3);
spct_n = size_x * size_y;
img_2d = reshape(img,[size_x*size_y,lambda_n]);
coeff = pca(img_2d);
% plot_coeff(coeff,5);
%%
%plot_pcs(img_2d,coeff,5,size_x,size_y) 
img_stack = zeros(size_x,size_y,clst_n);
for ii = 1:clst_n
    img_temp = img_2d*coeff(:,ii);
    img_temp = img_0_255(img_temp);
    img_stack(:,:,ii) = reshape(img_temp,size_x,size_y);
%     figure,imshow(reshape(img_temp,size_x,size_y))
%     title(['PC',num2str(ii)])
end
%%
[~,img_clsted] = max(img_stack,[],3);
img_color = psudo_clr(img_clsted,clst_n);           
img_color_re = reshape(img_color,[size_x,size_y,3]);
imshow(img_color_re)
end

function img_color = psudo_clr(img_clsted,clst_n)
    color_code = [0,0,0;0,0,1;0,1,1;0,1,0;1,0,0];
    img_color = zeros([size(img_clsted(:),1),3]);
    for kk = 1:clst_n
        locs = find(img_clsted(:) == kk);
        img_color(locs,:) = repmat(color_code(kk,:),length(locs),1);%color_code(kk,:);
    end
end

function plot_pcs(img_2d,coeff,num_PCs,size_x,size_y)
    for ii = 1:num_PCs
        img_temp = img_2d*coeff(:,ii);
        img_temp = img_0_255(img_temp);
        figure,imshow(reshape(img_temp,size_x,size_y))
        title(['PC',num2str(ii)])
    end
end

function img_new = img_0_255(img_temp)
    pix_0p005 = prctile(img_temp,0.5);
    pix_0p995 = prctile(img_temp,99.5);
    img_new = img_temp;
    img_new(img_new<pix_0p005) = pix_0p005;
    img_new(img_new>pix_0p995) = pix_0p995;
    img_new = (img_new- pix_0p005)/pix_0p995*255;
%     img_new = max(img_temp,pix_0p05);
%     img_new = min(img_new,pix_0p95);
end
    
function img_stack = read_tiff(img_name)
    info = imfinfo(img_name);
    num_imgs = length(info);
    img_stack = [];
    for kk = 1:num_imgs
        currentimg = imread(img_name,kk,'Info',info);
        img_stack(:,:,kk) = currentimg;
    end
end


function plot_coeff(coeff,num_PCs)
    for ii = 1:num_PCs
        figure,plot(coeff(:,ii))
        title(['PC',num2str(ii)])
    end
end
