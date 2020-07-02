%% gmm phasor function 
%% Reference 
% K.-C. Huang, J. Li, C. Zhang, and J.-X. Cheng, "SRS image cytometry for high-content single cell analysis," 
% in Multiphoton Microscopy in the Biomedical Sciences XIX (International Society for Optics and Photonics, 2019), Vol. 10882, p. 108822E.
% D. Fu and X. S. Xie, "Reliable Cell Segmentation Based on Spectral Phasor Analysis of Hyperspectral Stimulated Raman Scattering Imaging Data," 
% Anal. Chem. 86, 4115–4119 (2014).

%% initialize with manual phasor label
function gmm_label_3d = gmm_phasor_v1(SRS_re,label_phasor,display)
%% calculate phasor from srs stack
[size_x,size_y,nspec] = size(SRS_re);
fft_img = fft(permute(SRS_re,[3 1 2]));
fft_k1 = fft_img(2,:,:);
temp = fft_k1(:);
phasor = angle(temp);
phasor(:,2) = abs(temp);
phasor_xy = real(temp);
phasor_xy(:,2) = imag(temp);
phasor_xy(:,3) = 1:size(phasor_xy,1);
%% manual_phasor_label 
locs_cyto_man = find(sum(label_phasor,3)~=0);
% srs_flt_1 = mean(SRS_re,3).*reshape(locs_cyto_man,[size_x,size_y]);
%% get gmm initial
[gmm_ini,locs_all] = get_ini_phasor(phasor_xy,label_phasor);
PComponents = [length(locs_all{1}),length(locs_all{2}),length(locs_all{3}),length(locs_all{4})]/...
    sum([length(locs_all{1}),length(locs_all{2}),length(locs_all{3}),length(locs_all{4})]);
CLUSTER_NUM = 4;
if PComponents(1) == 0 %%% just in case where there is no lipid in the initial point -> how to fixxxxxxx??
    gmm_ini.Mu = gmm_ini.Mu(2:end,:);
    gmm_ini.Sigma = gmm_ini.Sigma(:,:,2:end);
    CLUSTER_NUM = 3;
    PComponents = PComponents(2:end);
end
S = struct('mu',gmm_ini.Mu,'Sigma',gmm_ini.Sigma,'ComponentProportion',PComponents);
options = statset('Display','final');

%% GMM model and plot
GMModel = fitgmdist([phasor_xy(locs_cyto_man,1),phasor_xy(locs_cyto_man,2)],CLUSTER_NUM,'Options',options,'Start',S);
if display == 1
    figure%,subplot(121)
    scatter(phasor_xy(locs_cyto_man,1),phasor_xy(locs_cyto_man,2))
    h = gca;
    hold on
    plot(GMModel.mu(:,1),GMModel.mu(:,2),'ok')
    fcontour(@(x1,x2)pdf(GMModel,[x1 x2]),[h.XLim h.YLim])%,'MeshDensity',100)
    color_ = {'r','g','y','b'};
    for ii = 1: CLUSTER_NUM
        x1 = GMModel.mu(ii,1) - 5:0.2:GMModel.mu(ii,1) + 5;
        x2 = GMModel.mu(ii,2) - 5:0.2:GMModel.mu(ii,2) + 5;
        [X1,X2] = meshgrid(x1,x2); X = [X1(:),X2(:)];
        y = mvnpdf(X,GMModel.mu(ii,:),abs(GMModel.Sigma(:,:,ii)));
        y = reshape(y,length(x2),length(x1));
        contour(x1,x2,y),hold on
    end
    title('phasor fitted with GMM-2')
end
%% get locs of all labels
locs_cyto_ = locs_cyto_man;%find(phasor(:,2) > threshold);
for nn = 1:CLUSTER_NUM
    mean_temp = GMModel.mu(nn,:);
    div_temp = sqrt(det(GMModel.Sigma(:,:,nn))/4/pi^2);
    inv_temp = inv(GMModel.Sigma(:,:,1));
    p_temp = [];
    for ii = 1:length(locs_cyto_)
        ind = locs_cyto_(ii);
        p_temp(ii) = exp(-0.5.*(phasor_xy(ind,1:2)-mean_temp)*inv_temp*(phasor_xy(ind,1:2)-mean_temp)')/div_temp;
    end
    if nn == 1
        p_all = p_temp';
    else
        p_all(:,nn) = p_temp';
    end
end

[~,indic] = max(p_all,[],2);

% figure,subplot(121)
label_mat = zeros(size_x*size_y,CLUSTER_NUM);
for nn = 1:CLUSTER_NUM
    locs_temp = locs_cyto_(find(indic == nn));
%     scatter(phasor_xy(locs_temp,1),phasor_xy(locs_temp,2),color_{nn}),hold on
    label_mat(locs_temp,nn) = 1;
end
% axis square,title('gmm cluster')
%% if no lipid label in the beginning -> no lipid label generated -> add zeros layer 
gmm_label_3d = reshape(label_mat,[size_x,size_y,CLUSTER_NUM]);
if size(gmm_label_3d,3) == 3
    gmm_label_3d = cat(3,zeros(size_x,size_y,1),gmm_label_3d);
% img_color_gmm = psudo_clr(label_mat,CLUSTER_NUM);
% subplot(122),imshow(img_color_gmm),title('gmm classification')
end

function [gmm_ini, locs_all] = get_ini_phasor(phasor_xy,label_phasor)
Mu = [];
locs_all = {};
Sigma = [];
for nn = 1:size(label_phasor,3) %% #class
    locs_temp = find(label_phasor(:,:,nn) == 1);%,[size(phasor_xy,1),1]);
    locs_all{nn} = locs_temp;
    if nn == 1
        Mu = mean(phasor_xy(locs_temp,1:2));
        Sigma = reshape(cov(phasor_xy(locs_temp,1:2)),[1,2,2]);
    else 
        Mu(nn,:) = mean(phasor_xy(locs_temp,1:2));
        Sigma(nn,:,:) = cov(phasor_xy(locs_temp,1:2));
    end
end
gmm_ini.Mu = Mu;
gmm_ini.Sigma = permute(Sigma,[2 3 1]);
end
end
    






