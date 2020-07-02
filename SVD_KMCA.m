%% k means cluster analysis 
%% Refenence
% M. Miljkovi?, T. Chernenko, M. J. Romeo, B. Bird, C. Matthäus, and M. Diem, "Label-free imaging of human cells: algorithms for image 
% reconstruction of Raman hyperspectral datasets," Analyst 135, 2002 (2010).

% clear,clc,close all
% clst_n = 5; %bkg, cytoplasma, nucleus, ER, lipid
% img_name ='001.tif';%read the hyperspectral stack 
% img = read_tiff(img_name);
% img = img(1:200,1:200,:);
function img_clsted = SVD_KMCA(img,clst_n)
size_x = size(img,1);
size_y = size(img,2);
lambda_n = size(img,3);
spct_n = size_x * size_y;
img_2d = reshape(img,[size_x*size_y,lambda_n]);
%% SVD
coeff = pca(img_2d);
img_2d_flt = zeros(size_x*size_y,clst_n);%clst_n
for ii = 1:clst_n
    img_temp = img_2d*coeff(:,ii);
    img_2d_flt(:,ii) = img_temp;
%     figure,imshow(reshape(img_temp,size_x,size_y))
%     title(['PC',num2str(ii)])
end
% random initialization
seed_1 = ceil(1000*rand(1,clst_n));
seed_clst = img_2d_flt(seed_1,:);
img_clsted = zeros(spct_n,1);
for iter = 1:500
    for ii = 1:spct_n
        dd = zeros(clst_n,1);
        for kk = 1:clst_n
            if ii == kk && iter == 1
                continue
            end
            dd(kk) = dist(img_2d_flt(ii,:),seed_clst(kk,:));
        end
        [y,loc] = min(dd);
        img_clsted(ii) = loc;
    end
    for kk = 1:clst_n
        locs_temp = find(img_clsted==kk);
        mean_temp = mean(img_2d_flt(locs_temp,:),1);
        seed_clst(kk,:) = mean_temp;
    end
end
img_clsted = reshape(img_clsted,[size_x,size_y]);
%%
img_color = psudo_clr(img_clsted,clst_n);
            
img_color_re = reshape(img_color,[size_x,size_y,3]);
imshow(img_color_re)
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


function dd = dist(spct1,spct_ctr)
    % for other distance function, refer to pdist 
    dd = sum((spct1-spct_ctr).^2);
end
%%
% function img_color = psudo_clr(img_clsted,clst_n)
%     color_code = [0,0,0;0,0,1;0,1,1;0,1,0;1,0,0];
%     img_color = zeros([size(img_clsted(:),1),3]);
%     for kk = 1:clst_n
%         locs = find(img_clsted(:) == kk);
%         img_color(locs,:) = repmat(color_code(kk,:),length(locs),1);%color_code(kk,:);
%     end
% end



