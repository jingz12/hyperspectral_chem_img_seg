%% reference:
% [1] P. Ghamisi, J. A. Benediktsson, and M. O. Ulfarsson, "Spectral–Spatial Classification of Hyperspectral Images Based on Hidden Markov Random Fields," 
%     IEEE Transactions on Geoscience and Remote Sensing 52, 2565–2574 (2014).
% [2] J. Li, J. M. Bioucas-Dias, and A. Plaza, "Spectral–Spatial Hyperspectral Image Segmentation Using Subspace Multinomial Logistic Regression and Markov Random Fields," 
%     IEEE Transactions on Geoscience and Remote Sensing 50, 809–823 (2012).
% [3] D. Fu and X. S. Xie, "Reliable Cell Segmentation Based on Spectral Phasor Analysis of Hyperspectral Stimulated Raman Scattering Imaging Data," 
%     Anal. Chem. 86, 4115–4119 (2014).

%% load hyperspectrum data and initial label
load('');
init_S = [lipid_spectrum,ER_spectrum,nuclei_spectrum,cyto_spectrum];% 2d, lambda-labels,
[~,nlabel] = size(init_S);
hspec_stack = exp_er2_srs_;%% 3d, X-Y-lambda
[size_x, size_y, nspec] = size(hspec_stack);

label_name = {'lipid','ER','nuclei','cytop'}; % load initial label
for ii = 1:4
    label_temp = imread([img_dir,'label\',label_name{ii},'.png']); %% each label is binary 
    label_temp(label_temp>0) = 255;
    if ii == 1
        label_crop = label_temp;
    else
        label_crop = cat(3,label_crop,label_temp);
    end
end
%% initialization
label_mat_ = [];label_spec_spat_mat_ = [];label_spec_spat_mat_nu_flt_ = [];
avg_mat_ = [];weighted_avg_mat_ = [];order_list = [];

kk_all_img = 0;
med_dim = 2; 
% parameter for segmentation 
alpha = [1,1,1,1];beta = [1, 1, 1, 1]; gama = [1,1,1,1];

SRS_re = permute(reshape(txt_temp,[size_x,nspec,size_y]),[1 3 2]);
SRS_2d = reshape(SRS_re,[size_x*size_y,nspec]);
[coeff2,score2] = pca(SRS_2d);
SRS_pca = reshape(score2(:,1:5),[size_x,size_y,5]);
%% filter the initial label if necessary 
label_crop_ = label_crop;
label_crop_(:,:,3) = medfilt2(label_crop(:,:,3),[med_dim med_dim]);
label_crop_(:,:,4) = logical(label_crop(:,:,4) + label_crop(:,:,3) - medfilt2(label_crop(:,:,3),[med_dim med_dim]));
%% spec spat label
label_phasor_2d = spec_spat_phasor(SRS_pca,label_crop,4,1,alpha,beta,gama);
label_phasor_2d_nuclei_flt = spec_spat_phasor(SRS_pca,label_crop_,4,1,alpha,beta,gama);
img_color_spec_spat = psudo_clr(label_phasor_2d,4);
%                 figure,imagesc(img_color_spec_spat),axis square,axis off
%%               
figure,subplot(231),imagesc(mean(SRS_re,3)),title('SRS unweighted mean'),axis square,axis off
subplot(233),imagesc(psudo_clr(label_crop,4)),title('phsor label'),axis square,axis off
subplot(234),imagesc(img_color_spec_spat),title('spec-spat-phasor label'),axis square,axis off
subplot(235),imagesc(label_crop_(:,:,3)),title('filtered nuclei label'),axis square,axis off
subplot(236),imagesc(psudo_clr(label_phasor_2d_nuclei_flt,4)),title(sprintf(...
    'spec-spat-phasor using filtered nuclei label,\nfilter size =%d',med_dim)),axis square,axis off
print(gcf,'-dpng',[saveDir,'\',num2str(nnn),'.png'])
    %% plot spec spat label with different parameters
%                 img_color_spec_spat_1_0p1 = psudo_clr(spec_spat_phasor(SRS_re,label_crop,10,1,0.1),4);
%                 img_color_spec_spat_0p1_0p1 = psudo_clr(spec_spat_phasor(SRS_re,label_crop,10,0.1,0.1),4);
%                 img_color_spec_spat_0p1_1 = psudo_clr(spec_spat_phasor(SRS_re,label_crop,10,0.1,1),4);
%                 img_color_spec_spat_0p01_0p1 = psudo_clr(spec_spat_phasor(SRS_re,label_crop,10,0.01,0.1),4);
    %% gmm phasor label
%                 label_gmm_phasor = gmm_phasor_v1(SRS_re,label_crop,1);
%                 img_color_gmm_phasor = psudo_clr(label_gmm_phasor,4);
%                 %%
%                 figure, subplot(221),imagesc(avg_crop),set(gca,'DataAspectRatio', [1 1 1]);
%                 axis off;title('average SRS')
%                 subplot(222),imagesc(img_color),set(gca,'DataAspectRatio', [1 1 1]);
%                 axis off;title('manual phasor label')
%                 subplot(223),imagesc(img_color_spec_spat),set(gca,'DataAspectRatio', [1 1 1]);
%                 axis off;title('spec+spat label')
%                 subplot(224),imagesc(img_color_gmm_phasor),set(gca,'DataAspectRatio', [1 1 1]);
%                 axis off;title('GMM phasor label')
%                 print(gcf,'-dpng',[myDir,'label\merge_crop_2\',num2str(nnn),'.png'])
%                 save([myDir,'label\raw_mat\',num2str(nnn),'.mat'],'SRS_re')

%% save parameters 
para_name = {'alpha','beta','gama'};
txt_name = [saveDir,'\','para.txt'];
for ii = 1:length(para_name)
    fid = fopen(txt_name,'a');
    fprintf(fid,[para_name{ii},'\n'])
    fclose(fid);
    eval(sprintf(['data_temp = ',para_name{ii},';']));
    dlmwrite(txt_name,data_temp,'-append','delimiter','\t','precision',6);
end

