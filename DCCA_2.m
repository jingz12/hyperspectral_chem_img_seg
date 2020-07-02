%% divisive correlation cluster analysis DCCA
%% Refenence
% M. Miljkovi?, T. Chernenko, M. J. Romeo, B. Bird, C. Matthäus, and M. Diem, "Label-free imaging of human cells: algorithms for image 
% reconstruction of Raman hyperspectral datasets," Analyst 135, 2002 (2010).
% the current version have to chose the divisive cluster manually
clear,clc,close all
clst_n = 5; %bkg, cytoplasma, nucleus, ER, lipid
img_name ='001.tif';%read the hyperspectral stack 
img = read_tiff(img_name);
img = img(1:20,1:20,:);
size_x = size(img,1);
size_y = size(img,2);
lambda_n = size(img,3);
spct_n = size_x * size_y;
global img_2d
img_2d = reshape(img,[size_x*size_y,lambda_n]);
%%
%img_2d = rand(10,4);
dist_1 = pdist(img_2d,'correlation');
dist_1_sqr = squareform(dist_1);
[x,y] = find(dist_1_sqr == max(dist_1));
if length(x) > 1
    x = x(1);y = y(1);
end
S_L = img_2d(x,:);
S_M = img_2d(y,:);
%S_other = img_2d;S_other(x) = [];S_other(y) = [];
n_iter = 1;
%[locs1,locs2] = divisive(S_L,S_M,S_other);
locs = 1:spct_n;
%clst = zeros(spct_n,1);
global clsted nn doRun
clsted = {};
nn = 1;
doRun = true;
divisive(S_L,S_M,img_2d,locs);

%% 
function divisive(S_L,S_M,img_2d_,locs)
    global clsted nn img_2d doRun
    if nn > 6
        doRun = false;
        return%error('Breaking out of function divisive');%return
    end
    corr_L_temp = corr(img_2d_',S_L');
    corr_M_temp = corr(img_2d_',S_M');
    locs_M = corr_L_temp < corr_M_temp;
    locs_L = 1 - locs_M;
    locs_M_temp = locs.*locs_M'; locs_M_temp(locs_M_temp==0)=[];
    locs_L_temp = locs.*locs_L'; locs_L_temp(locs_L_temp==0)=[];
%     clst(locs_M_temp) = n_iter*2 - 1;
%     clst(locs_L_temp) = n_iter*2;
    clsted{2*nn} = locs_M_temp;
    clsted{2*nn-1} = locs_L_temp;
    nn = nn + 1;    
    parfor iii = 1:2
        if  doRun == true
            if iii == 1
                locs_temp = locs.*locs_M';locs_temp(locs_temp==0) = [];
                S_temp = img_2d(locs_temp,:);
            else
                locs_temp = locs.*locs_L';locs_temp(locs_temp==0) = [];
                S_temp = img_2d(locs_temp,:);
            end
            dist_temp = pdist(S_temp,'correlation');
            dist_1_sqr = squareform(dist_temp);
            [x,y] = find(dist_1_sqr == max(dist_temp));
            if length(x) >= 2
                x = x(1);y = y(1);
            end
    %         S_L = S_temp(x,:);
    %         S_M = S_temp(y,:);
            divisive(S_temp(x,:),S_temp(y,:),S_temp,locs_temp);
            %iii = iii - 1;
        end
    end
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


%% 
% clear, clc 
% z = fib(1,2)
% global var nn 
% var = [];
% nn = 1;
% function z = fib(x,y)
%     global var nn
%     z = x+y;
%     if z >100
%         return
%     end
%     if nn == 1 
%         var = z;
%     else
%         var(nn) = z;
%     end   
%     z = fib(y,z)% call with diff arguments, or the loop is infinited 
%     nn = nn + 1;
% end


