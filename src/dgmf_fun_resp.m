function dfx = dgmf_fun_resp(x)
% Gather global variables form the problem initialization
global K delta beta
global rows cols old_rows old_cols offset_x offset_y
global fft_img training_set_size training_set

global gmf_resps
global GMF_resp
global GMF_exps
global Sum_GMF_exps
global alphas_ks

    % Unpack the GMF parameters from the x vector:
    sigma = x(1);
    L = x(2);
    T = x(3);
    w = x(4);
    b = x(5);

    %% Compute the GMF response:
    U = linspace(-cols/2, cols/2-1, cols);
    V = linspace(-rows/2, rows/2-1, rows);

    [u, v] = meshgrid(U, V);

    dgmf_resps_padded.dsigma = zeros(rows, cols, K);
    dgmf_resps_padded.dL = zeros(rows, cols, K);
    dgmf_resps_padded.dT = zeros(rows, cols, K);

    % Acummulate the gradients using the chain rule:
    dfx = zeros(training_set_size*old_rows*old_cols, 5);
    
    %% Precompute the GMF kernels:
    fft_gmf.dsigma = cell(K,1);
    fft_gmf.dL = cell(K,1);
    fft_gmf.dT = cell(K,1);
    dgmf = gmf_gradients(u, v, sigma, L, T, K, delta);

    for k = 1:K
        fft_gmf.dsigma{k} = fft2(dgmf.dsigma{k});
        fft_gmf.dL{k} = fft2(dgmf.dL{k});
        fft_gmf.dT{k} = fft2(dgmf.dT{k});
    end

    %% Apply to all the images from the training set
    img_idx = 1:(old_rows*old_cols);
    for img_i = 1:training_set_size
        %% Apply the GMF filter to the current image:        
        for k=1:K
            dgmf_resps_padded.dsigma(:,:,k) = fftshift(ifft2(fft_gmf.dsigma{k} .* fft_img(:,:,training_set(img_i))));
            dgmf_resps_padded.dL(:,:,k) = fftshift(ifft2(fft_gmf.dL{k} .* fft_img(:,:,training_set(img_i))));
            dgmf_resps_padded.dT(:,:,k) = fftshift(ifft2(fft_gmf.dT{k} .* fft_img(:,:,training_set(img_i))));
        end

        % Trim the responses to fit in the original image dimensions
        dGMF_resp_k_dsigma = dgmf_resps_padded.dsigma(offset_y:(offset_y+old_rows-1),offset_x:(offset_x+old_cols-1),:);
        dGMF_resp_k_dL = dgmf_resps_padded.dL(offset_y:(offset_y+old_rows-1),offset_x:(offset_x+old_cols-1),:);
        dGMF_resp_k_dT = dgmf_resps_padded.dT(offset_y:(offset_y+old_rows-1),offset_x:(offset_x+old_cols-1),:);
        
        %% Compute the gradients with respect to the alphas computation
        dGMF_exps_dsigma = beta*GMF_exps(:,:,:,img_i).*dGMF_resp_k_dsigma;
        dGMF_exps_dL = beta*GMF_exps(:,:,:,img_i).*dGMF_resp_k_dL;
        dGMF_exps_dT = beta*GMF_exps(:,:,:,img_i).*dGMF_resp_k_dT;
        
        dSum_GMF_exps_dsigma = sum(dGMF_exps_dsigma, 3);
        dSum_GMF_exps_dL = sum(dGMF_exps_dL, 3);
        dSum_GMF_exps_dT = sum(dGMF_exps_dT, 3);
        
        Sum_GMF_exps_2 = Sum_GMF_exps(:,:,img_i).^2;
        
        dalphas_ks_dsigma = bsxfun(@rdivide, bsxfun(@times, Sum_GMF_exps(:,:,img_i), dGMF_exps_dsigma) - bsxfun(@times, GMF_exps(:,:,:,img_i), dSum_GMF_exps_dsigma), Sum_GMF_exps_2);
        dalphas_ks_dL = bsxfun(@rdivide, bsxfun(@times, Sum_GMF_exps(:,:,img_i), dGMF_exps_dL) - bsxfun(@times, GMF_exps(:,:,:,img_i), dSum_GMF_exps_dL), Sum_GMF_exps_2);
        dalphas_ks_dT = bsxfun(@rdivide, bsxfun(@times, Sum_GMF_exps(:,:,img_i), dGMF_exps_dT) - bsxfun(@times, GMF_exps(:,:,:,img_i), dSum_GMF_exps_dT), Sum_GMF_exps_2);
        
        dfx(img_idx, 1) = reshape(sum(alphas_ks(:,:,:,img_i).*dGMF_resp_k_dsigma + gmf_resps(:,:,:,img_i).*dalphas_ks_dsigma, 3), [old_rows*old_cols, 1]);
        dfx(img_idx, 2) = reshape(sum(alphas_ks(:,:,:,img_i).*dGMF_resp_k_dL + gmf_resps(:,:,:,img_i).*dalphas_ks_dL, 3), [old_rows*old_cols, 1]);
        dfx(img_idx, 3) = reshape(sum(alphas_ks(:,:,:,img_i).*dGMF_resp_k_dT + gmf_resps(:,:,:,img_i).*dalphas_ks_dT, 3), [old_rows*old_cols, 1]);

        img_idx = img_idx + old_rows*old_cols;
    end
    dfx(:,4) = GMF_resp(:);
    dfx(:,5) = 1;
end