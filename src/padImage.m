function padded_img = padImage(img)
% Pad the original image to reach a 2 powered size. The new space generated
% is filled with symmetric view of the image.
%
%	Fernando Cervantes Sanchez, June 2019
%	iie.fercer@gmail.com
    
    [rows, cols] = size(img);
    nearest_2p_dim = ceil(log2(rows));
    rows_2 = 2^nearest_2p_dim;
    cols_2 = 2^nearest_2p_dim;
    
    padded_img = padarray(img, [rows_2/2 - rows/2, cols_2/2-cols/2], 'symmetric');
end
