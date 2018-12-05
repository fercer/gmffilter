/***********************************************************************************************************************************
CENTRO DE INVESTIGACION EN MATEMATICAS
DOCTORADO EN CIENCIAS DE LA COMPUTACION
FERNANDO CERVANTES SANCHEZ

FILE NAME : gmf_extension.cpp

PURPOSE : Implements the functions for the GMF filter method implemented in c++ for its use within python.

FILE REFERENCES :
Name        I / O        Description
None----       ----------

ABNORMAL TERMINATION CONDITIONS, ERROR AND WARNING MESSAGES :
None
************************************************************************************************************************************/

#include "gmf.h"


/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv C implementations: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

/************************************************************************************************************
* FUNCTION NAME: filterGMF_impl
*
* PURPOSE: Computes the GMF filtering of the img_src using the parameters passed. 
* This function uses the fast convolution in the Fourier space of frequencies.
*
* ARGUMENTS:
* ARGUMENT                  TYPE                      I/O  DESCRIPTION
* img_src                   double *                   I   Source image to be filtered
* img_dst                   double *                   O   Response of the source image to the GMF
* ang_dst                   double *                   O   Angles of max response of the image to the GMF
* height                    const int                  I   Height of the source image
* width                     const int                  I   Width of the source image
* par_T                     const int                  I   Gaussian curve trails cut point (for the GMF template)
* par_L                     const int                  I   Height of the GMF template
* par_sigma                 const double               I   Spread of the Gaussian curve (for the GMF template)
* par_K                     const int                  I   Orientations to rotate and apply the template
*
* RETURNS:
* The response of the GMF filter applied to the source image, and the angles of max response at each pixel.
*
************************************************************************************************************/
void filterGMF_impl(double* img_src, double* img_dst, double* ang_dst, const int height, const int width, const int par_T, const int par_L, const double par_sigma, const int par_K)
{	
	/* Compute the nearest 2 powered dimension where the input image fits */
	const unsigned int max_dim = (height > width) ? height : width;
	const double remainder_2p_dim = log2(max_dim) - floor(log2(max_dim));
	const unsigned int nearest_2p_dim = (unsigned int)pow(2, floor(log2(max_dim)) + (remainder_2p_dim > 1e-12 ? 1 : 0));
	
	/* Offset for the zero padded image */
	const unsigned int GMF_kernel_height = par_L + 6;
	const unsigned int GMF_kernel_width = par_T + 3;
	const unsigned int offset_y = GMF_kernel_height / 2;
	const unsigned int offset_x = GMF_kernel_width / 2;

	ft_complex * fft_zp_img;
	ft_complex ** fft_GMF_kernels;
	ft_complex * fft_resp_to_GMF_kernel;

	double* resp_to_GMF_kernels_data;
	double* max_resp;
	double* max_resp_angle;

	/* ----------------------------------------------------------- Memory allocation ----------------------------------------------------------- */
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	/* Zero pad the input images: */
	fft_zp_img = allocate_ft_complex((nearest_2p_dim / 2 + 1)*nearest_2p_dim);

	double* zp_img = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
	double* img_data_ptr = img_src;
	for (unsigned int i = 0; i < height; i++)
	{
		for (unsigned int j = 0; j < width; j++, img_data_ptr++)
		{			
			*(zp_img + i * nearest_2p_dim + j) = *img_data_ptr;
		}
	}

	ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_zp_img);
	ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_zp_img);
	ft_release_forward;

	free(zp_img);	
	

	fft_GMF_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_GMF_kernels + k) = allocate_ft_complex((nearest_2p_dim / 2 + 1) * nearest_2p_dim);
	}

	/* Arrays for filtering process */
	resp_to_GMF_kernels_data = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	max_resp = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	max_resp_angle = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	fft_resp_to_GMF_kernel = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim / 2 + 1));

	/* ----------------------------------------------------------- Filters definition ----------------------------------------------------------- */
	generateGMFTemplate_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim);
	
	/* ----------------------------------------------------------- Apply the GMF filter --------------------------------------------------------- */
	double minima_filter_response = MY_INF;
	double maxima_filter_response = -MY_INF;
	
	/* First Kernel: */
	for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
	{
		const double real_result = ft_real(*(fft_GMF_kernels)+i) * ft_real(fft_zp_img + i) -
			ft_imag(*(fft_GMF_kernels)+i) * ft_imag(fft_zp_img + i);
		const double imag_result = ft_real(*(fft_GMF_kernels)+i) * ft_imag(fft_zp_img + i) +
			ft_real(fft_zp_img + i) * ft_imag(*(fft_GMF_kernels)+i);

		ft_real_assign(fft_resp_to_GMF_kernel + i) = real_result;
		ft_imag_assign(fft_resp_to_GMF_kernel + i) = imag_result;

		*(max_resp + i) = -MY_INF;
		*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = -MY_INF;
	}

	ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_resp_to_GMF_kernel, resp_to_GMF_kernels_data);
	ft_backward(nearest_2p_dim, nearest_2p_dim, fft_resp_to_GMF_kernel, resp_to_GMF_kernels_data);
	ft_release_backward;
		
	for (unsigned int k = 1; k < par_K; k++)
	{
		for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
		{
			const double real_result = ft_real(*(fft_GMF_kernels + k) + i) * ft_real(fft_zp_img + i) -
				ft_imag(*(fft_GMF_kernels + k) + i) * ft_imag(fft_zp_img + i);

			const double imag_result = ft_real(*(fft_GMF_kernels + k) + i) * ft_imag(fft_zp_img + i) + ft_real(fft_zp_img + i) * ft_imag(*(fft_GMF_kernels + k) + i);

			ft_real_assign(fft_resp_to_GMF_kernel + i) = real_result;
			ft_imag_assign(fft_resp_to_GMF_kernel + i) = imag_result;

			/* Check if the response to the previous kernel is higher to the current response: */
			if (*(max_resp + i) < *(resp_to_GMF_kernels_data + i))
			{
				*(max_resp + i) = *(resp_to_GMF_kernels_data + i);
				*(max_resp_angle + i) = MY_PI * (double)(k - 1) / (double)par_K;
			}

			if (*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) < *(resp_to_GMF_kernels_data + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim))
			{
				*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = *(resp_to_GMF_kernels_data + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim);
				*(max_resp_angle + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = MY_PI * (double)(k - 1) / (double)par_K;
			}
		}

		ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_resp_to_GMF_kernel, resp_to_GMF_kernels_data);
		ft_backward(nearest_2p_dim, nearest_2p_dim, fft_resp_to_GMF_kernel, resp_to_GMF_kernels_data);
		ft_release_backward;
		
	}

	for (unsigned int i = 0; i < nearest_2p_dim*nearest_2p_dim; i++)
	{
		/* Check if the response to the last kernel is higher to the current response: */
		if (*(max_resp + i) < *(resp_to_GMF_kernels_data + i))
		{
			*(max_resp + i) = *(resp_to_GMF_kernels_data + i);
			*(max_resp_angle + i) = MY_PI * (double)(par_K - 1) / (double)par_K;
		}
		
		*(max_resp + i) = *(max_resp + i) / (double)(nearest_2p_dim*nearest_2p_dim);
	}


	double* GMF_resp_it = img_dst;
	double* GMF_resp_angles_it = ang_dst;
	for (unsigned int i = 0; i < height; i++)
	{
		for (unsigned int j = 0; j < width; j++, GMF_resp_it++, GMF_resp_angles_it++)
		{
			*GMF_resp_it = *(max_resp + (i + offset_y) * nearest_2p_dim + offset_x + j);
			*GMF_resp_angles_it = *(max_resp_angle + (i + offset_y) * nearest_2p_dim + offset_x + j);
		}
	}
	
	/* ----------------------------------------------------------- Memory deallocation ----------------------------------------------------------- */
	deallocate_ft_complex(fft_zp_img);

	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(fft_GMF_kernels + k));
	}
	deallocate_ft_complex(fft_GMF_kernels);

	free(resp_to_GMF_kernels_data);
	free(max_resp);
	free(max_resp_angle);

	deallocate_ft_complex(fft_resp_to_GMF_kernel);
	ft_close;
}



void generateGMFTemplate_impl(ft_complex ** fft_GMF_kernels, const int par_T, const int par_L, const double par_sigma, const int par_K, const int nearest_2p_dim)
{	
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	
	const unsigned int GMF_kernel_height = par_L + 6;
	const unsigned int GMF_kernel_width = par_T + 3;

	/* Generate the GMF kernels */
	double *GMF_base_kernel_temp = (double*)calloc(GMF_kernel_width, sizeof(double));
	double Gauss_val, curr_x = -floor((double)(par_T) / 2.0);
	double sum_GMF_kernel_base = 0.0;

	double *GMF_kernel_it = GMF_base_kernel_temp + 1;
	for (unsigned int j = 0; j < (par_T + 1 - par_T % 2); j++, GMF_kernel_it++, curr_x += 1.0)
	{
		Gauss_val = 1.0 - exp(-curr_x * curr_x / (2.0 * par_sigma * par_sigma));
		*GMF_kernel_it = Gauss_val;
		sum_GMF_kernel_base += Gauss_val;
	}

	const double mean_GMF_kernel_base = sum_GMF_kernel_base / (double)par_T;
	GMF_kernel_it = GMF_base_kernel_temp + 1;
	for (unsigned int j = 0; j < par_T; j++, GMF_kernel_it++)
	{
		*GMF_kernel_it = (*GMF_kernel_it - mean_GMF_kernel_base) / (sum_GMF_kernel_base * par_L);
	}

	/* Orientate the K GMF kernels */
	/* Transform the zero padded kernels to the frequencies space */
	double *GMF_base_kernel = (double*)calloc(GMF_kernel_height * GMF_kernel_width, sizeof(double));
	/* At 0 degrees: */
	for (unsigned int i = 0; i < par_L; i++)
	{
		memcpy(GMF_base_kernel + (i + 3) * GMF_kernel_width + 1, GMF_base_kernel_temp + 1, par_T * sizeof(double));
	}
	free(GMF_base_kernel_temp);

	double *zp_GMF_rotated_kernel = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
	double *GMF_rotated_kernel = (double*)malloc(GMF_kernel_height * GMF_kernel_width * sizeof(double));
	
	for (unsigned int k = 0; k < par_K; k++)
	{
		rotateBicubic_impl(GMF_base_kernel, GMF_rotated_kernel, (double)k * 180.0 / (double)par_K, GMF_kernel_height, GMF_kernel_width);
		
		for (unsigned int i = 0; i < GMF_kernel_height; i++)
		{
			memcpy(zp_GMF_rotated_kernel + i * nearest_2p_dim, GMF_rotated_kernel + i*GMF_kernel_width, GMF_kernel_width * sizeof(double));
		}
						
		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_GMF_rotated_kernel, *(fft_GMF_kernels + k));
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_GMF_rotated_kernel, *(fft_GMF_kernels + k));
		ft_release_forward;
	}
	

	free(zp_GMF_rotated_kernel);
	free(GMF_rotated_kernel);
	free(GMF_base_kernel);
	ft_close;
}



/************************************************************************************************************
* FUNCTION NAME: filterGMF_impl
*
* PURPOSE: Computes the GMF filtering of the img_src using the parameters passed.
* This function uses the fast convolution in the Fourier space of frequencies.
*
* ARGUMENTS:
* ARGUMENT                  TYPE                      I/O  DESCRIPTION
* img_src                   double *                   I   Source image to be filtered
* img_dst                   double *                   O   Response of the source image to the GMF
* ang_dst                   double *                   O   Angles of max response of the image to the GMF
* height                    const int                  I   Height of the source image
* width                     const int                  I   Width of the source image
* par_T                     const int                  I   Gaussian curve trails cut point (for the GMF template)
* par_L                     const int                  I   Height of the GMF template
* par_sigma                 const double               I   Spread of the Gaussian curve (for the GMF template)
* par_K                     const int                  I   Orientations to rotate and apply the template
*
* RETURNS:
* The response of the GMF filter applied to the source image, and the angles of max response at each pixel.
*
************************************************************************************************************/
void filterGMFUntrimmed_impl(double* img_src, double* img_dst, double* ang_dst, const int height, const int width, const int par_T, const int par_L, const double par_sigma, const int par_K)
{
	/* Compute the nearest 2 powered dimension where the input image fits */
	const unsigned int max_dim = (height > width) ? height : width;
	const double remainder_2p_dim = log2(max_dim) - floor(log2(max_dim));
	const unsigned int nearest_2p_dim = (unsigned int)pow(2, floor(log2(max_dim)) + (remainder_2p_dim > 1e-12 ? 1 : 0));

	/* Offset for the zero padded image */
	const unsigned int GMF_kernel_height = par_L + 6;
	const unsigned int GMF_kernel_width = par_T + 3;
	const unsigned int GMF_kernel_size = (int)ceil(sqrt(GMF_kernel_height*GMF_kernel_height + GMF_kernel_width*GMF_kernel_width) / 2.0);
	const unsigned int offset = GMF_kernel_size / 2;

	ft_complex * fft_zp_img;
	ft_complex ** fft_GMF_kernels;
	ft_complex * fft_resp_to_GMF_kernel;

	double *resp_to_GMF_kernels_data;
	double *max_resp;
	double *max_resp_angle;

	/* ----------------------------------------------------------- Memory allocation ----------------------------------------------------------- */
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	/* Zero pad the input images: */
	fft_zp_img = allocate_ft_complex((nearest_2p_dim / 2 + 1)*nearest_2p_dim);

	double* zp_img = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
	double* img_data_ptr = img_src;
	for (unsigned int i = 0; i < height; i++)
	{
		for (unsigned int j = 0; j < width; j++, img_data_ptr++)
		{
			*(zp_img + i * nearest_2p_dim + j) = *(double*)img_data_ptr;
		}
	}

	ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_zp_img);
	ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_zp_img);
	ft_release_forward;
	
	free(zp_img);

	fft_GMF_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_GMF_kernels + k) = allocate_ft_complex((nearest_2p_dim / 2 + 1) * nearest_2p_dim);
	}

	/* Arrays for filtering process */
	resp_to_GMF_kernels_data = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	max_resp = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	max_resp_angle = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	fft_resp_to_GMF_kernel = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim / 2 + 1));

	/* ----------------------------------------------------------- Filters definition ----------------------------------------------------------- */
	generateGMFTemplate_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim);

	/* ----------------------------------------------------------- Apply the GMF filter --------------------------------------------------------- */
	double minima_filter_response = MY_INF;
	double maxima_filter_response = -MY_INF;

	/* First Kernel: */
	for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
	{
		const double real_result = ft_real(*(fft_GMF_kernels)+i) * ft_real(fft_zp_img + i) -
			ft_imag(*(fft_GMF_kernels)+i) * ft_imag(fft_zp_img + i);

		const double imag_result = ft_real(*(fft_GMF_kernels)+i) * ft_imag(fft_zp_img + i) + ft_real(fft_zp_img + i) * ft_imag(*(fft_GMF_kernels)+i);

		ft_real_assign(fft_resp_to_GMF_kernel + i) = real_result;
		ft_imag_assign(fft_resp_to_GMF_kernel + i) = imag_result;

		*(max_resp + i) = -MY_INF;
		*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = -MY_INF;
	}

	ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_resp_to_GMF_kernel, resp_to_GMF_kernels_data);
	ft_backward(nearest_2p_dim, nearest_2p_dim, fft_resp_to_GMF_kernel, resp_to_GMF_kernels_data);
	ft_release_backward;
	
	for (unsigned int k = 1; k < par_K; k++)
	{
		for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
		{
			const double real_result = ft_real(*(fft_GMF_kernels + k) + i) * ft_real(fft_zp_img + i) -
				ft_imag(*(fft_GMF_kernels + k) + i) * ft_imag(fft_zp_img + i);

			const double imag_result = ft_real(*(fft_GMF_kernels + k) + i) * ft_imag(fft_zp_img + i) + ft_real(fft_zp_img + i) * ft_imag(*(fft_GMF_kernels + k) + i);

			ft_real_assign(fft_resp_to_GMF_kernel + i) = real_result;
			ft_real_assign(fft_resp_to_GMF_kernel + i) = imag_result;

			/* Check if the response to the previous kernel is higher to the current response: */
			if (*(max_resp + i) < *(resp_to_GMF_kernels_data + i))
			{
				*(max_resp + i) = *(resp_to_GMF_kernels_data + i);
				*(max_resp_angle + i) = MY_PI * (double)(k - 1) / (double)par_K;
			}

			if (*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) < *(resp_to_GMF_kernels_data + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim))
			{
				*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = *(resp_to_GMF_kernels_data + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim);
				*(max_resp_angle + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = MY_PI * (double)(k - 1) / (double)par_K;
			}
		}

		ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_resp_to_GMF_kernel, resp_to_GMF_kernels_data);
		ft_backward(nearest_2p_dim, nearest_2p_dim, fft_resp_to_GMF_kernel, resp_to_GMF_kernels_data);
		ft_release_backward;
	}

	for (unsigned int i = 0; i < nearest_2p_dim*nearest_2p_dim; i++)
	{
		/* Check if the response to the last kernel is higher to the current response: */
		if (*(max_resp + i) < *(resp_to_GMF_kernels_data + i))
		{
			*(max_resp + i) = *(resp_to_GMF_kernels_data + i);
			*(max_resp_angle + i) = MY_PI * (double)(par_K - 1) / (double)par_K;
		}

		*(max_resp + i) = *(max_resp + i) / (double)(nearest_2p_dim*nearest_2p_dim);
	}

	double* GMF_resp_it = img_dst;
	double* GMF_resp_angles_it = ang_dst;
	for (unsigned int i = 0; i < height; i++)
	{
		for (unsigned int j = 0; j < width; j++, GMF_resp_it++, GMF_resp_angles_it++)
		{
			*GMF_resp_it = *(max_resp + (i + offset) * nearest_2p_dim + offset + j);
			*GMF_resp_angles_it = *(max_resp_angle + (i + offset) * nearest_2p_dim + offset + j);
		}
	}

	/* ----------------------------------------------------------- Memory deallocation ----------------------------------------------------------- */
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(fft_GMF_kernels + k));
	}
	deallocate_ft_complex(fft_GMF_kernels);

	free(resp_to_GMF_kernels_data);
	free(max_resp);
	free(max_resp_angle);

	deallocate_ft_complex(fft_resp_to_GMF_kernel);
	ft_close;
}




void generateGMFTemplateUntrimmed_impl(ft_complex ** fft_GMF_kernels, const int par_T, const int par_L, const double par_sigma, const int par_K, const int nearest_2p_dim)
{
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	
	const unsigned int GMF_kernel_height = par_L + 6;
	const unsigned int GMF_kernel_width = par_T + 3;
	const unsigned int GMF_kernel_size = (int)ceil(sqrt(GMF_kernel_height * GMF_kernel_height + GMF_kernel_width * GMF_kernel_width) / 2.0);
	const unsigned int GMF_kernel_offset_x = GMF_kernel_size - GMF_kernel_width / 2;
	const unsigned int GMF_kernel_offset_y = GMF_kernel_size - GMF_kernel_height / 2;

	/* Generate the GMF kernels */
	double* GMF_base_kernel_temp = (double*)calloc(GMF_kernel_width, sizeof(double));
	double Gauss_val, curr_x = -floor((double)(par_T) / 2.0);
	double sum_GMF_kernel_base = 0.0;

	double* GMF_kernel_it = GMF_base_kernel_temp + 1;
	for (unsigned int j = 0; j < (par_T + 1 - par_T % 2); j++, GMF_kernel_it++, curr_x += 1.0)
	{
		Gauss_val = 1.0 - exp(-curr_x * curr_x / (2.0 * par_sigma * par_sigma));
		*GMF_kernel_it = Gauss_val;
		sum_GMF_kernel_base += Gauss_val;
	}

	const double mean_GMF_kernel_base = sum_GMF_kernel_base / (double)par_T;
	GMF_kernel_it = GMF_base_kernel_temp + 1;
	for (unsigned int j = 0; j < par_T; j++, GMF_kernel_it++)
	{
		*GMF_kernel_it = (*GMF_kernel_it - mean_GMF_kernel_base) / (sum_GMF_kernel_base * par_L);
	}

	/* Orientate the K GMF kernels */
	/* Transform the zero padded kernels to the frequencies space */
	double *GMF_base_kernel = (double*)calloc(GMF_kernel_size * GMF_kernel_size, sizeof(double));

	/* At 0 degrees: */
	for (unsigned int i = 0; i < par_L; i++)
	{
		memcpy(GMF_base_kernel + (i + 3 + GMF_kernel_offset_y) * GMF_kernel_width + 1 + GMF_kernel_offset_x, GMF_base_kernel_temp + 1, par_T * sizeof(double));
	}
	free(GMF_base_kernel_temp);

	double *zp_GMF_rotated_kernel = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
	double *GMF_rotated_kernel = (double*)calloc(GMF_kernel_height * GMF_kernel_width, sizeof(double));
	for (unsigned int k = 0; k < par_K; k++)
	{
		rotateBicubic_impl(GMF_base_kernel, GMF_rotated_kernel, (double)k * 180.0 / (double)par_K, GMF_kernel_height, GMF_kernel_width);

		for (unsigned int i = 0; i < GMF_kernel_height; i++)
		{
			memcpy(zp_GMF_rotated_kernel + i * nearest_2p_dim, GMF_rotated_kernel + i*GMF_kernel_width, GMF_kernel_width * sizeof(double));
		}

		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_GMF_rotated_kernel, *(fft_GMF_kernels + k));
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_GMF_rotated_kernel, *(fft_GMF_kernels + k));
		ft_release_forward;
	}

	free(zp_GMF_rotated_kernel);
	free(GMF_rotated_kernel);
	free(GMF_base_kernel);
	ft_close;
}




/************************************************************************************************************
* FUNCTION NAME: bicubicInterpolation
*
* PURPOSE: Calculates the interpolation between points in 'src' using bicubic interponalation of the known points
*
* ARGUMENTS:
* ARGUMENT                  TYPE                      I/O  DESCRIPTION
* src                       double*                    I   Source matrix of points
* x                         const double               I   Position to interpolate in the x axis
* y                         const double               I   Position to interpolate in the y axis
*
* RETURNS:
* The interpolated value in the (x, y) position.
*
************************************************************************************************************/
inline double bicubicInterpolation_impl(double* src, const double x, const double y, const unsigned int src_height, const unsigned int src_width)
{
	const int x_int = (int)floor(x);
	const int y_int = (int)floor(y);

	const double s0_x = x - (double)x_int;
	const double s_1_x = 1.0 + s0_x;
	const double s1_x = 1.0 - s0_x;
	const double s2_x = 2.0 - s0_x;

	const double s0_y = y - (double)y_int;
	const double s_1_y = 1.0 + s0_y;
	const double s1_y = 1.0 - s0_y;
	const double s2_y = 2.0 - s0_y;

	/* Compute coefficients for the x-axis interpolation */
	const double ux_1 = -0.5 * s_1_x*s_1_x*s_1_x + 2.5 * s_1_x * s_1_x - 4.0 * s_1_x + 2.0;
	const double ux0 = 1.5 * s0_x*s0_x*s0_x - 2.5 * s0_x*s0_x + 1.0;
	const double ux1 = 1.5 * s1_x*s1_x*s1_x - 2.5 * s1_x*s1_x + 1.0;
	const double ux2 = -0.5 * s2_x*s2_x*s2_x + 2.5 * s2_x * s2_x - 4.0 * s2_x + 2.0;

	/* Compute coefficients for the y-axis interpolation */
	const double uy_1 = -0.5 * s_1_y*s_1_y*s_1_y + 2.5 * s_1_y * s_1_y - 4.0 * s_1_y + 2.0;
	const double uy0 = 1.5 * s0_y*s0_y*s0_y - 2.5 * s0_y*s0_y + 1.0;
	const double uy1 = 1.5 * s1_y*s1_y*s1_y - 2.5 * s1_y*s1_y + 1.0;
	const double uy2 = -0.5 * s2_y*s2_y*s2_y + 2.5 * s2_y * s2_y - 4.0 * s2_y + 2.0;

	return
		(*(src + (y_int - 1) * src_width + x_int - 1) * uy_1 +
			*(src + y_int * src_width + x_int - 1) * uy0 +
			*(src + (y_int + 1) * src_width + x_int - 1) * uy1 +
			*(src + (y_int + 2) * src_width + x_int - 1) * uy2) * ux_1 +
			(*(src + (y_int - 1) * src_width + x_int) * uy_1 +
				*(src + y_int * src_width + x_int) * uy0 +
				*(src + (y_int + 1) * src_width + x_int) * uy1 +
				*(src + (y_int + 2) * src_width + x_int) * uy2) * ux0 +
				(*(src + (y_int - 1) * src_width + x_int + 1) * uy_1 +
					*(src + y_int * src_width + x_int + 1) * uy0 +
					*(src + (y_int + 1) * src_width + x_int + 1) * uy1 +
					*(src + (y_int + 2) * src_width + x_int + 1) * uy2) * ux1 +
					(*(src + (y_int - 1) * src_width + x_int + 2) * uy_1 +
						*(src + y_int * src_width + x_int + 2) * uy0 +
						*(src + (y_int + 1) * src_width + x_int + 2) * uy1 +
						*(src + (y_int + 2) * src_width + x_int + 2) * uy2) * ux2;
}





/************************************************************************************************************
* FUNCTION NAME: rotateBicubic
*
* PURPOSE: Rotates the matrix 'src' by 'theta' radian degrees, using the bicubic interpolation.
*
* ARGUMENTS:
* ARGUMENT                  TYPE                      I/O  DESCRIPTION
* src                       double*                    I   Source matrix of points
* dst                       double*                    O   Destination matrix
* theta                     const double               I   Degrees the source matrix is rotated
* src_height                const unsigned int         I   Source height
* src_width                 const unsigned int         I   Source width
*
* RETURNS:
* The matrix 'src' rotated 'theta' degrees into the 'dst' matrix.
*
************************************************************************************************************/
void rotateBicubic_impl(double* src, double* dst, const double theta, const int src_height, const int src_width)
{
	double c_theta, s_theta;
	double half_src_height, half_src_width;

	if (theta == 0.0)
	{
		memcpy(dst, src, src_width * src_height * sizeof(double));
		return;
	}
	else if (theta == 90.0)
	{
		c_theta = 0.0;
		s_theta = 1.0;
		half_src_height = floor((double)src_height / 2.0);
		half_src_width = floor((double)src_width / 2.0);
	}
	else
	{
		c_theta = cos(theta / 180.0 * MY_PI);
		s_theta = sin(theta / 180.0 * MY_PI);

		half_src_height = (double)(src_height - 1) / 2.0;
		half_src_width = (double)(src_width - 1) / 2.0;
	}

	/* Copy the source image to a temporal zero padded matrix. */
	const int off_set = (int)floor(sqrt(half_src_height * half_src_height / 4.0 + half_src_width*half_src_width / 4.0) + 4.0);
	double *src_temp = (double*)calloc((src_height + 2 * off_set) * (src_width + 2 * off_set), sizeof(double));
	for (int i = 0; i < src_height; i++)
	{
		memcpy(src_temp + (src_width + 2 * off_set) * (i + off_set) + off_set, src + i*src_width, src_width * sizeof(double));
	}
	
	for (int i = 0; i < src_height; i++)
	{
		for (int j = 0; j < src_width; j++)
		{		
			const double src_x = c_theta * ((double)j - half_src_width) - s_theta * ((double)i - half_src_height) + half_src_width;
			const double src_y = s_theta * ((double)j - half_src_width) + c_theta * ((double)i - half_src_height) + half_src_height;
			if (src_x < -(double)off_set || src_x >= (double)(src_width + off_set) || src_y < -(double)off_set || src_y >= (double)(src_height + off_set))
			{
				*(dst + src_width*i + j) = 0.0;
			}
			else
			{
				*(dst + src_width*i + j) = bicubicInterpolation_impl(src_temp, src_x + (double)off_set, src_y + (double)off_set, src_height + 2 * off_set, src_width + 2 * off_set);
			}
		}
	}

	free(src_temp);
}

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ C++ implementations: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/



/*vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Python interface: vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv*/

/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Python interface: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^*/
