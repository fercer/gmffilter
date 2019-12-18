function fx = gmf_fun_resp(x)
% Gather global variables form the problem initialization
global K delta beta
global rows cols old_rows old_cols offset_x offset_y
global fft_img training_set_size training_set

global gmf_resps
global GMF_exps
global Sum_GMF_exps
global alphas_ks
global GMF_resp

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

    gmf_resps = zeros(old_rows, old_cols, K, training_set_size);
    GMF_resp = zeros(old_rows, old_cols, training_set_size);
    GMF_exps = zeros(old_rows, old_cols, K, training_set_size);
    Sum_GMF_exps = zeros(old_rows, old_cols, K, training_set_size);
    alphas_ks = zeros(old_rows, old_cols, K, training_set_size);
    gmf_resps_padded = zeros(rows, cols, K);

    %% Precompute the GMF kernels:
    fft_gmf = cell(K,1);
    gmf = gmf_functions(u, v, sigma, L, T, K, delta);
    for k = 1:K
        fft_gmf{k} = fft2(gmf{k});
    end

    %% Apply to all the images from the training set
    for img_i = 1:training_set_size
        %% Apply the GMF filter to the current image:        
        for k=1:K
            gmf_resps_padded(:,:,k) = fftshift(ifft2(fft_gmf{k} .* fft_img(:,:,training_set(img_i))));
        end

        % Trim the respnses to fit in the original image dimensions
        gmf_resps(:,:,:,img_i) = gmf_resps_padded(offset_y:(offset_y+old_rows-1),offset_x:(offset_x+old_cols-1),:);

        %% Compute the weights of each orientation as the probability of observing the vessel pixel with the corresponding orientation
        GMF_exps(:,:,:,img_i) = exp(beta*gmf_resps(:,:,:,img_i));
        Sum_GMF_exps(:,:,img_i) = sum(GMF_exps(:,:,:,img_i), 3);
        alphas_ks(:,:,:,img_i) = bsxfun(@rdivide, GMF_exps(:,:,:,img_i), Sum_GMF_exps(:,:,img_i));
        GMF_resp(:, :, img_i) = sum(alphas_ks(:,:,:,img_i).*gmf_resps(:,:,:,img_i),3);
    end
    
    fx = w*GMF_resp(:) + b;
end