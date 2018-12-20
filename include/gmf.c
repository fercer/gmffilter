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


double bicubicInterpolation_impl(double* src, const double x, const double y, const unsigned int src_height, const unsigned int src_width)
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


void generateGMFTemplate_impl(ft_complex ** fft_GMF_kernels, const int par_T, const int par_L, const double par_sigma, const int par_K, const int nearest_2p_dim, unsigned int *GMF_kernel_height, unsigned int *GMF_kernel_width)
{
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	
	*GMF_kernel_height = par_L + 6;
	*GMF_kernel_width = par_T + 3;

	/* Generate the GMF kernels */
	double *GMF_base_kernel_temp = (double*)calloc(*GMF_kernel_width, sizeof(double));
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
	double *GMF_base_kernel = (double*)calloc(*GMF_kernel_height * *GMF_kernel_width, sizeof(double));
	/* At 0 degrees: */
	for (unsigned int i = 0; i < par_L; i++)
	{
		memcpy(GMF_base_kernel + (i + 3) * *GMF_kernel_width + 1, GMF_base_kernel_temp + 1, par_T * sizeof(double));
	}
	free(GMF_base_kernel_temp);

	double *zp_GMF_rotated_kernel = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
	double *GMF_rotated_kernel = (double*)malloc(*GMF_kernel_height * *GMF_kernel_width * sizeof(double));
	
	for (unsigned int k = 0; k < par_K; k++)
	{
		rotateBicubic_impl(GMF_base_kernel, GMF_rotated_kernel, (double)k * 180.0 / (double)par_K, *GMF_kernel_height, *GMF_kernel_width);
		
		for (unsigned int i = 0; i < *GMF_kernel_height; i++)
		{
			memcpy(zp_GMF_rotated_kernel + i * nearest_2p_dim, GMF_rotated_kernel + i* *GMF_kernel_width, *GMF_kernel_width * sizeof(double));
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


void generateGMFTemplateUntrimmed_impl(ft_complex ** fft_GMF_kernels, const int par_T, const int par_L, const double par_sigma, const int par_K, const int nearest_2p_dim, unsigned int * GMF_kernel_height, unsigned int * GMF_kernel_width)
{
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	
	*GMF_kernel_height = par_L + 6;
	*GMF_kernel_width = par_T + 3;
	const unsigned int GMF_kernel_size = (int)ceil(sqrt(*GMF_kernel_height * *GMF_kernel_height + *GMF_kernel_width * *GMF_kernel_width) / 2.0);
	const unsigned int GMF_kernel_offset_x = GMF_kernel_size - *GMF_kernel_width / 2;
	const unsigned int GMF_kernel_offset_y = GMF_kernel_size - *GMF_kernel_height / 2;

	/* Generate the GMF kernels */
	double* GMF_base_kernel_temp = (double*)calloc(*GMF_kernel_width, sizeof(double));
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
		memcpy(GMF_base_kernel + (i + 3 + GMF_kernel_offset_y) * *GMF_kernel_width + 1 + GMF_kernel_offset_x, GMF_base_kernel_temp + 1, par_T * sizeof(double));
	}
	free(GMF_base_kernel_temp);

	double *zp_GMF_rotated_kernel = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
	double *GMF_rotated_kernel = (double*)calloc(*GMF_kernel_height * *GMF_kernel_width, sizeof(double));
	for (unsigned int k = 0; k < par_K; k++)
	{
		rotateBicubic_impl(GMF_base_kernel, GMF_rotated_kernel, (double)k * 180.0 / (double)par_K, *GMF_kernel_height, *GMF_kernel_width);

		for (unsigned int i = 0; i < *GMF_kernel_height; i++)
		{
			memcpy(zp_GMF_rotated_kernel + i * nearest_2p_dim, GMF_rotated_kernel + i* *GMF_kernel_width, *GMF_kernel_width * sizeof(double));
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


void applyGMFWithAngles(ft_complex * fft_img_src, double* img_dst, double* ang_dst, const int nearest_2p_dim, const int height, const int width, const int GMF_kernel_height, const int GMF_kernel_width, const int par_K, ft_complex ** fft_GMF_kernels)
{	
	/* Offset for the zero padded image */
	const unsigned int offset_y = GMF_kernel_height / 2;
	const unsigned int offset_x = GMF_kernel_width / 2;

	ft_variables(nearest_2p_dim, nearest_2p_dim);

	/* Arrays for filtering process */
	double * resp_to_GMF_kernels_data = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	double * max_resp = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	double * max_resp_angle = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	ft_complex * fft_convolution_GMF_kernel = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim / 2 + 1));

	/* First Kernel: */
	for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
	{
		const double real_result = ft_real(*(fft_GMF_kernels)+i) * ft_real(fft_img_src + i) -
			ft_imag(*(fft_GMF_kernels)+i) * ft_imag(fft_img_src + i);
		const double imag_result = ft_real(*(fft_GMF_kernels)+i) * ft_imag(fft_img_src + i) +
			ft_real(fft_img_src + i) * ft_imag(*(fft_GMF_kernels)+i);

		ft_real_assign(fft_convolution_GMF_kernel + i) = real_result;
		ft_imag_assign(fft_convolution_GMF_kernel + i) = imag_result;

		*(max_resp + i) = -MY_INF;
		*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = -MY_INF;
	}

	ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
	ft_backward(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
	ft_release_backward;
		
	for (unsigned int k = 1; k < par_K; k++)
	{
		for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
		{
			const double real_result = ft_real(*(fft_GMF_kernels + k) + i) * ft_real(fft_img_src + i) -
				ft_imag(*(fft_GMF_kernels + k) + i) * ft_imag(fft_img_src + i);

			const double imag_result = ft_real(*(fft_GMF_kernels + k) + i) * ft_imag(fft_img_src + i) + 
				ft_real(fft_img_src + i) * ft_imag(*(fft_GMF_kernels + k) + i);

			ft_real_assign(fft_convolution_GMF_kernel + i) = real_result;
			ft_imag_assign(fft_convolution_GMF_kernel + i) = imag_result;

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

		ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
		ft_backward(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
		ft_release_backward;
		
	}
	deallocate_ft_complex(fft_convolution_GMF_kernel);
	
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
	free(resp_to_GMF_kernels_data);
	ft_close;

	for (unsigned int i = 0; i < height; i++)
	{
		memcpy(img_dst + i*width, max_resp + (i + offset_y) * nearest_2p_dim + offset_x, width*sizeof(double));
		memcpy(ang_dst + i*width, max_resp_angle + (i + offset_y) * nearest_2p_dim + offset_x, width*sizeof(double));
	}
	
	free(max_resp);
	free(max_resp_angle);
}


void applyGMF(ft_complex * fft_img_src, double* img_dst, const int nearest_2p_dim, const int height, const int width, const int GMF_kernel_height, const int GMF_kernel_width, const int par_K, ft_complex ** fft_GMF_kernels)
{	
	/* Offset for the zero padded image */
	const unsigned int offset_y = GMF_kernel_height / 2;
	const unsigned int offset_x = GMF_kernel_width / 2;

	ft_variables(nearest_2p_dim, nearest_2p_dim);

	/* Arrays for filtering process */
	double * resp_to_GMF_kernels_data = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	double * max_resp = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	ft_complex * fft_convolution_GMF_kernel = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim / 2 + 1));

	/* First Kernel: */
	for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
	{
		const double real_result = ft_real(*(fft_GMF_kernels)+i) * ft_real(fft_img_src + i) -
			ft_imag(*(fft_GMF_kernels)+i) * ft_imag(fft_img_src + i);
		const double imag_result = ft_real(*(fft_GMF_kernels)+i) * ft_imag(fft_img_src + i) +
			ft_real(fft_img_src + i) * ft_imag(*(fft_GMF_kernels)+i);

		ft_real_assign(fft_convolution_GMF_kernel + i) = real_result;
		ft_imag_assign(fft_convolution_GMF_kernel + i) = imag_result;

		*(max_resp + i) = -MY_INF;
		*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = -MY_INF;
	}

	ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
	ft_backward(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
	ft_release_backward;
		
	for (unsigned int k = 1; k < par_K; k++)
	{
		for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
		{
			const double real_result = ft_real(*(fft_GMF_kernels + k) + i) * ft_real(fft_img_src + i) -
				ft_imag(*(fft_GMF_kernels + k) + i) * ft_imag(fft_img_src + i);

			const double imag_result = ft_real(*(fft_GMF_kernels + k) + i) * ft_imag(fft_img_src + i) + 
				ft_real(fft_img_src + i) * ft_imag(*(fft_GMF_kernels + k) + i);

			ft_real_assign(fft_convolution_GMF_kernel + i) = real_result;
			ft_imag_assign(fft_convolution_GMF_kernel + i) = imag_result;

			/* Check if the response to the previous kernel is higher to the current response: */
			if (*(max_resp + i) < *(resp_to_GMF_kernels_data + i))
			{
				*(max_resp + i) = *(resp_to_GMF_kernels_data + i);
			}

			if (*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) < *(resp_to_GMF_kernels_data + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim))
			{
				*(max_resp + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim) = *(resp_to_GMF_kernels_data + i + (nearest_2p_dim / 2 - 1) * nearest_2p_dim);
			}
		}

		ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
		ft_backward(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
		ft_release_backward;
		
	}
	deallocate_ft_complex(fft_convolution_GMF_kernel);
	
	for (unsigned int i = 0; i < nearest_2p_dim*nearest_2p_dim; i++)
	{
		/* Check if the response to the last kernel is higher to the current response: */
		if (*(max_resp + i) < *(resp_to_GMF_kernels_data + i))
		{
			*(max_resp + i) = *(resp_to_GMF_kernels_data + i);
		}
		
		*(max_resp + i) = *(max_resp + i) / (double)(nearest_2p_dim*nearest_2p_dim);
	}
	free(resp_to_GMF_kernels_data);
	ft_close;

	for (unsigned int i = 0; i < height; i++)
	{
		memcpy(img_dst + i*width, max_resp + (i + offset_y) * nearest_2p_dim + offset_x, width*sizeof(double));
	}
	
	free(max_resp);
}

void singleScaleGMFilter(double * raw_input, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, const double par_sigma, const unsigned int par_K, const unsigned char untrimmed_kernels)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	for (unsigned int y = 0; y < height; y++)
	{
        memcpy(zp_img + y*nearest_2p_dim, raw_input + y*width, width * sizeof(double));
	}
		
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	ft_complex* fft_img_src = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
	
	/* Perform the Fourier Transform: */
	ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
	ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
	ft_release_forward;
	
	free(zp_img);
	
	// Compute the filter kernels:
	ft_complex** fft_GMF_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_GMF_kernels+k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2+1));
	}
	
	unsigned int GMF_kernel_height, GMF_kernel_width;
	if (untrimmed_kernels)
	{
		generateGMFTemplateUntrimmed_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
	}
	else
	{
		generateGMFTemplate_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
	}
	
	// Apply the single-scale filter:
	applyGMF(fft_img_src, output, nearest_2p_dim, height, width, GMF_kernel_height, GMF_kernel_width, par_K, fft_GMF_kernels);

	// Apply the mask, the mask(x, y) must be {0, 1}:
	char * m_ptr = mask;
	double *o_ptr = output;
        
	for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
	{
		*o_ptr = *o_ptr * (double)*m_ptr;
	}
	
	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(fft_GMF_kernels + k));
	}
	deallocate_ft_complex(fft_GMF_kernels);
	deallocate_ft_complex(fft_img_src);	
	ft_close;
}


void singleScaleGMFilterWithAngles(double * raw_input, char * mask, double * output, double * angles_output, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, const double par_sigma, const unsigned int par_K, const unsigned char untrimmed_kernels)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	for (unsigned int y = 0; y < height; y++)
	{
        memcpy(zp_img + y*nearest_2p_dim, raw_input + y*width, width * sizeof(double));
    }
    
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	ft_complex* fft_img_src = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
	
	/* Perform the Fourier Transform: */
	ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
	ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
	ft_release_forward;
	
	free(zp_img);
	
	// Compute the filter kernels:
	ft_complex** fft_GMF_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_GMF_kernels+k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2+1));
	}
	
	unsigned int GMF_kernel_height, GMF_kernel_width;
	if (untrimmed_kernels)
	{
		generateGMFTemplateUntrimmed_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
	}
	else
	{
		generateGMFTemplate_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
	}
	
	// Apply the single-scale filter:
	applyGMFWithAngles(fft_img_src, output, angles_output, nearest_2p_dim, height, width, GMF_kernel_height, GMF_kernel_width, par_K, fft_GMF_kernels);

	// Apply the mask, the mask(x, y) must be {0, 1}:
	char * m_ptr = mask;
	double *o_ptr = output;
	double *a_ptr = angles_output;
    
	for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++, a_ptr++)
	{
		*o_ptr = *o_ptr * (double)*m_ptr;
		*a_ptr = *a_ptr * (double)*m_ptr;
	}
	
	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(fft_GMF_kernels + k));
	}
	deallocate_ft_complex(fft_GMF_kernels);
	deallocate_ft_complex(fft_img_src);	
	ft_close;
}


void singleScaleGMFilter_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, const double par_sigma, const unsigned int par_K, const unsigned char untrimmed_kernels)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	ft_complex* fft_img_src = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
	
	// Compute the filter kernels:
	ft_complex** fft_GMF_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_GMF_kernels+k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2+1));
	}
	
	unsigned int GMF_kernel_height, GMF_kernel_width;
	if (untrimmed_kernels)
	{
		generateGMFTemplateUntrimmed_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
	}
	else
	{
		generateGMFTemplate_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
	}
	
	// Apply the single-scale filter:
	for (unsigned int i = 0; i < n_inputs; i++)
	{
		for (unsigned int y = 0; y < height; y++)
		{
            memcpy(zp_img + y*nearest_2p_dim, raw_input + i*width*height + y*width, width * sizeof(double));
		}
			
		/* Perform the Fourier Transform: */
		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
		ft_release_forward;
		
		applyGMF(fft_img_src, output + i*height*width, nearest_2p_dim, height, width, GMF_kernel_height, GMF_kernel_width, par_K, fft_GMF_kernels);

		// Apply the mask, the mask(x, y) must be {0, 1}:
		char * m_ptr = mask + i*height*width;
		double *o_ptr = output + i*height*width;
		
		for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
		{
			*o_ptr = *o_ptr * (double)*m_ptr;
		}
	}
	
	free(zp_img);
	
	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(fft_GMF_kernels + k));
	}
	deallocate_ft_complex(fft_GMF_kernels);
	deallocate_ft_complex(fft_img_src);	
	ft_close;
}


void singleScaleGMFilterWithAngles_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, double * angles_output, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, const double par_sigma, const unsigned int par_K, const unsigned char untrimmed_kernels)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	ft_complex* fft_img_src = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
	
	// Compute the filter kernels:
	ft_complex** fft_GMF_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_GMF_kernels+k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2+1));
	}
	
	unsigned int GMF_kernel_height, GMF_kernel_width;
	if (untrimmed_kernels)
	{
		generateGMFTemplateUntrimmed_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
	}
	else
	{
		generateGMFTemplate_impl(fft_GMF_kernels, par_T, par_L, par_sigma, par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
	}
	
	// Apply the single-scale filter:
	for (unsigned int i = 0; i < n_inputs; i++)
	{
		for (unsigned int y = 0; y < height; y++)
		{
            memcpy(zp_img + y*nearest_2p_dim, raw_input + i*width*height + y*width, width * sizeof(double));
		}
			
		/* Perform the Fourier Transform: */
		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
		ft_release_forward;
		
		applyGMFWithAngles(fft_img_src, output + i*height*width, angles_output + i*height*width, nearest_2p_dim, height, width, GMF_kernel_height, GMF_kernel_width, par_K, fft_GMF_kernels);

		// Apply the mask, the mask(x, y) must be {0, 1}:
		char * m_ptr = mask + i*height*width;
		double *o_ptr = output + i*height*width;
		double *a_ptr = angles_output + i*height*width;
		
		for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++, a_ptr++)
		{
			*o_ptr = *o_ptr * (double)*m_ptr;
			*a_ptr = *a_ptr * (double)*m_ptr;
		}
	}
	
	free(zp_img);
	
	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(fft_GMF_kernels + k));
	}
	deallocate_ft_complex(fft_GMF_kernels);
	deallocate_ft_complex(fft_img_src);	
	ft_close;
}


void multiscaleGMFilter(double * raw_input, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char untrimmed_kernels)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	for (unsigned int y = 0; y < height; y++)
	{
        memcpy(zp_img + y*nearest_2p_dim, raw_input + y*width, width * sizeof(double));
	}
	
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	ft_complex* fft_img_src = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
	
	/* Perform the Fourier Transform: */
	ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
	ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
	ft_release_forward;
	
	free(zp_img);
	
	// Compute the filter kernels:
	unsigned int GMF_kernel_height, GMF_kernel_width;
	ft_complex** fft_GMF_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_GMF_kernels + k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2+1));
	}
	
	for (unsigned int s = 0; s < sigma_scales; s++)
	{		
		if (untrimmed_kernels)
		{
			generateGMFTemplateUntrimmed_impl(fft_GMF_kernels, par_T, par_L, *(par_sigma + s), par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
		}
		else
		{
			generateGMFTemplate_impl(fft_GMF_kernels, par_T, par_L, *(par_sigma + s), par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
		}
	
		// Apply the single-scale filter:
		applyGMF(fft_img_src, output + s*height*width, nearest_2p_dim, height, width, GMF_kernel_height, GMF_kernel_width, par_K, fft_GMF_kernels);

		// Apply the mask, the mask(x, y) must be {0, 1}:
		char * m_ptr = mask;
		double *o_ptr = output;
		
		for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
		{
			*o_ptr = *o_ptr * (double)*m_ptr;
		}
	}
	
	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(fft_GMF_kernels + k));
	}
	deallocate_ft_complex(fft_GMF_kernels);
	deallocate_ft_complex(fft_img_src);	
	ft_close;
}


void multiscaleGMFilterWithAngles(double * raw_input, char * mask, double * output, double * angles_output, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char untrimmed_kernels)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	for (unsigned int y = 0; y < height; y++)
	{
        memcpy(zp_img + y*nearest_2p_dim, raw_input + y*width, width * sizeof(double));
	}
	
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	ft_complex* fft_img_src = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
	
	/* Perform the Fourier Transform: */
	ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
	ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
	ft_release_forward;
	
	free(zp_img);
	
	// Compute the filter kernels:
	unsigned int GMF_kernel_height, GMF_kernel_width;
	ft_complex** fft_GMF_kernels = allocate_ft_complex_pointers(par_K);
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_GMF_kernels + k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2+1));
	}
	
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		
		if (untrimmed_kernels)
		{
			generateGMFTemplateUntrimmed_impl(fft_GMF_kernels, par_T, par_L, *(par_sigma + s), par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
		}
		else
		{
			generateGMFTemplate_impl(fft_GMF_kernels, par_T, par_L, *(par_sigma + s), par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
		}
	
		// Apply the single-scale filter:
		applyGMFWithAngles(fft_img_src, output + s*height*width, angles_output + s*height*width, nearest_2p_dim, height, width, GMF_kernel_height, GMF_kernel_width, par_K, fft_GMF_kernels);

		// Apply the mask, the mask(x, y) must be {0, 1}:
		char * m_ptr = mask;
		double *o_ptr = output;
		double *a_ptr = angles_output;
		
		for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++, a_ptr++)
		{
			*o_ptr = *o_ptr * (double)*m_ptr;
			*a_ptr = *a_ptr * (double)*m_ptr;
		}
	}
	
	// Free all memory used for image the filtering:
	for (unsigned int k = 0; k < par_K; k++)
	{
		deallocate_ft_complex(*(fft_GMF_kernels + k));
	}
	deallocate_ft_complex(fft_GMF_kernels);
	deallocate_ft_complex(fft_img_src);	
	ft_close;
}


void multiscaleGMFilter_multipleinputs(double * raw_input, const unsigned int n_inputs, char * mask, double * output, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char untrimmed_kernels)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	ft_complex* fft_img_src = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
	
	// Compute the filter kernels:
	unsigned int GMF_kernel_height, GMF_kernel_width;
	ft_complex*** fft_GMF_kernels = allocate_ft_complex_pointers_of_pointers(sigma_scales);
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		*(fft_GMF_kernels + s) = allocate_ft_complex_pointers(par_K);
		for (unsigned int k = 0; k < par_K; k++)
		{
			*(*(fft_GMF_kernels+s)+k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2+1));
		}
	
		if (untrimmed_kernels)
		{
			generateGMFTemplateUntrimmed_impl(*(fft_GMF_kernels + s), par_T, par_L, *(par_sigma + s), par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
		}
		else
		{
			generateGMFTemplate_impl(*(fft_GMF_kernels + s), par_T, par_L, *(par_sigma + s), par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
		}
	}
	
	// Apply the single-scale filter:
	for (unsigned int i = 0; i < n_inputs; i++)
	{
		for (unsigned int y = 0; y < height; y++)
		{
            memcpy(zp_img + y*nearest_2p_dim, raw_input + i*width*height + y*width, width * sizeof(double));
		}
			
		/* Perform the Fourier Transform: */
		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
		ft_release_forward;
	
		for (unsigned int s = 0; s < sigma_scales; s++)
		{
			applyGMF(fft_img_src, output + s*height*width + i*sigma_scales*height*width, nearest_2p_dim, height, width, GMF_kernel_height, GMF_kernel_width, par_K, *(fft_GMF_kernels + s));

			// Apply the mask, the mask(x, y) must be {0, 1}:
			char * m_ptr = mask + i*height*width;
			double *o_ptr = output + s*height*width + i*sigma_scales*height*width;
		
			for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++)
			{
				*o_ptr = *o_ptr * (double)*m_ptr;
			}
		}
	}
	
	free(zp_img);
	
	// Free all memory used for image the filtering:
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		for (unsigned int k = 0; k < par_K; k++)
		{
			deallocate_ft_complex(*(*(fft_GMF_kernels + s) + k));
		}
		deallocate_ft_complex(*(fft_GMF_kernels + s));
	}
	deallocate_ft_complex(fft_GMF_kernels);
	deallocate_ft_complex(fft_img_src);	
	ft_close;
}


void multiscaleGMFilterWithAngles_multipleinputs(double * raw_input, unsigned int n_inputs, char * mask, double * output, double * angles_output, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char untrimmed_kernels)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	ft_complex* fft_img_src = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
	
	// Compute the filter kernels:
	unsigned int GMF_kernel_height, GMF_kernel_width;
	ft_complex*** fft_GMF_kernels = allocate_ft_complex_pointers_of_pointers(sigma_scales);
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		*(fft_GMF_kernels + s) = allocate_ft_complex_pointers(par_K);
		for (unsigned int k = 0; k < par_K; k++)
		{
			*(*(fft_GMF_kernels+s)+k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2+1));
		}
	
		if (untrimmed_kernels)
		{
			generateGMFTemplateUntrimmed_impl(*(fft_GMF_kernels + s), par_T, par_L, *(par_sigma + s), par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
		}
		else
		{
			generateGMFTemplate_impl(*(fft_GMF_kernels + s), par_T, par_L, *(par_sigma + s), par_K, nearest_2p_dim, &GMF_kernel_height, &GMF_kernel_width);
		}
	}
	
	// Apply the single-scale filter:
	for (unsigned int i = 0; i < n_inputs; i++)
	{
		for (unsigned int y = 0; y < height; y++)
		{
            memcpy(zp_img + y*nearest_2p_dim, raw_input + i*width*height + y*width, width * sizeof(double));
		}
			
		/* Perform the Fourier Transform: */
		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, fft_img_src);
		ft_release_forward;
	
		for (unsigned int s = 0; s < sigma_scales; s++)
		{
			applyGMFWithAngles(fft_img_src, output + s*height*width + i*sigma_scales*height*width, angles_output + s*height*width + i*sigma_scales*height*width, nearest_2p_dim, height, width, GMF_kernel_height, GMF_kernel_width, par_K, *(fft_GMF_kernels + s));

			// Apply the mask, the mask(x, y) must be {0, 1}:
			char * m_ptr = mask + i*height*width;
			double *o_ptr = output + s*height*width + i*sigma_scales*height*width;
			double *a_ptr = angles_output + s*height*width + i*sigma_scales*height*width;
		
			for (unsigned int xy = 0; xy < height*width; xy++, m_ptr++, o_ptr++, a_ptr++)
			{
				*o_ptr = *o_ptr * (double)*m_ptr;
				*a_ptr = *a_ptr * (double)*m_ptr;
			}
		}
	}
	
	free(zp_img);
	
	// Free all memory used for image the filtering:
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		for (unsigned int k = 0; k < par_K; k++)
		{
			deallocate_ft_complex(*(*(fft_GMF_kernels + s) + k));
		}
		deallocate_ft_complex(*(fft_GMF_kernels + s));
	}
	deallocate_ft_complex(fft_GMF_kernels);
	deallocate_ft_complex(fft_img_src);	
	ft_close;
}



#ifdef BUILDING_PYTHON_MODULE
static PyObject* gmfFilter(PyObject *self, PyObject *args)
{
    PyArrayObject *raw_input;
    char *raw_input_data = NULL;
    npy_intp raw_input_stride, n_imgs = 1;
    npy_intp height, width;
    
    PyArrayObject *multiscale_par_sigma = NULL;
    char * par_sigma_data = NULL;
    npy_intp par_sigma_stride, par_sigma_scales = 1;
    
    PyArrayObject *mask = NULL;
    char * mask_data = NULL;
    npy_intp mask_stride;
    
    unsigned int par_T;
    unsigned int par_L;
    unsigned int par_K;
    
    DEBMSG("Filtering by Gaussian matched filters");
    
    if (!PyArg_ParseTuple(args, "O!IIO!IO!", &PyArray_Type, &raw_input, &par_T, &par_L, &PyArray_Type, &multiscale_par_sigma, &par_K, &PyArray_Type, &mask))
    {
        return NULL;
    }
    
    par_sigma_data = ((PyArrayObject*)multiscale_par_sigma)->data;
    par_sigma_scales = ((PyArrayObject*)multiscale_par_sigma)->dimensions[0];
    par_sigma_stride = ((PyArrayObject*)multiscale_par_sigma)->strides[((PyArrayObject*)multiscale_par_sigma)->nd-1];
    
    DEBNUMMSG("T = %i", par_T);
    DEBNUMMSG(", L = %i", par_L);
    DEBNUMMSG(", K = %i", par_K);
    DEBNUMMSG(", n sigma scales = %i\n", (intpar_sigma_scales);
    
    if (((PyArrayObject*)raw_input)->nd > 2)
    {
        n_imgs = ((PyArrayObject*)raw_input)->dimensions[0];
        height = ((PyArrayObject*)raw_input)->dimensions[1];
        width = ((PyArrayObject*)raw_input)->dimensions[2];
    }
    else
    {
        n_imgs = 1;
        height = ((PyArrayObject*)raw_input)->dimensions[0];
        width = ((PyArrayObject*)raw_input)->dimensions[1];
    }
    
    raw_input_data  = ((PyArrayObject*)raw_input)->data;
    raw_input_stride = ((PyArrayObject*)raw_input)->strides[((PyArrayObject*)raw_input)->nd - 1];
    
	mask_data = ((PyArrayObject*)mask)->data;
	mask_stride = ((PyArrayObject*)mask)->strides[((PyArrayObject*)mask)->nd - 1];
 
    
    PyObject * gmf_response = NULL;
    if (n_imgs > 1 && par_sigma_scales > 1) 
    {
        npy_intp gmf_response_shape[] = { n_imgs, par_sigma_scales, height, width };      
        gmf_response = PyArray_SimpleNew(4, &gmf_response_shape[0], NPY_DOUBLE);
    }
    else if (n_imgs > 1 && par_sigma_scales == 1)
    {
        npy_intp gmf_response_shape[] = { n_imgs, height, width };        
        gmf_response = PyArray_SimpleNew(3, &gmf_response_shape[0], NPY_DOUBLE);
    }
    else if (n_imgs == 1 && par_sigma_scales > 1)
    {
        npy_intp gmf_response_shape[] = { par_sigma_scales, height, width };      
        gmf_response = PyArray_SimpleNew(3, &gmf_response_shape[0], NPY_DOUBLE);
    }
    else
    {
        npy_intp gmf_response_shape[] = { height, width };        
        gmf_response = PyArray_SimpleNew(2, &gmf_response_shape[0], NPY_DOUBLE);
    }
    
    char * gmf_response_data = ((PyArrayObject*)gmf_response)->data;
    npy_intp gmf_response_stride = ((PyArrayObject*)gmf_response)->strides[((PyArrayObject*)gmf_response)->nd-1];
    
    if (n_imgs > 1)
    {
        DEBMSG("Multiscale Gaussian matched filtering over multiple images\n");
        multiscaleGMFilter_multipleinputs((double*)raw_input_data, n_imgs, mask_data, (double*)gmf_response_data, height, width, par_T, par_L, (double*)par_sigma_data, par_sigma_scales, par_K, 0);
    }
    else
    {
        DEBMSG("Multiscale Gaussian matched filtering over a single image\n");
        multiscaleGMFilter((double*)raw_input_data, mask_data, (double*)gmf_response_data, height, width, par_T, par_L, (double*)par_sigma_data, par_sigma_scales, par_K, 0);
    }
      
    return gmf_response;
}
#endif



#ifdef BUILDING_PYTHON_MODULE
static PyObject* gmfFilterWithAngles(PyObject *self, PyObject *args)
{
    PyArrayObject *raw_input;
    char *raw_input_data = NULL;
    npy_intp raw_input_stride, n_imgs = 1;
    npy_intp height, width;
    
    PyArrayObject *multiscale_par_sigma = NULL;
    char * par_sigma_data = NULL;
    npy_intp par_sigma_stride, par_sigma_scales = 1;
    
    PyArrayObject *mask = NULL;
    char * mask_data = NULL;
    npy_intp mask_stride;
    
    unsigned int par_T;
    unsigned int par_L;
    unsigned int par_K;
    
    if (!PyArg_ParseTuple(args, "O!IIO!IO!", &PyArray_Type, &raw_input, &par_T, &par_L, &PyArray_Type, &multiscale_par_sigma, &par_K, &PyArray_Type, &mask))
    {
        return NULL;
    }
    
    par_sigma_data = ((PyArrayObject*)multiscale_par_sigma)->data;
    par_sigma_scales = ((PyArrayObject*)multiscale_par_sigma)->dimensions[0];
    par_sigma_stride = ((PyArrayObject*)multiscale_par_sigma)->strides[((PyArrayObject*)multiscale_par_sigma)->nd-1];
        
    if (((PyArrayObject*)raw_input)->nd > 2)
    {
        n_imgs = ((PyArrayObject*)raw_input)->dimensions[0];
        height = ((PyArrayObject*)raw_input)->dimensions[1];
        width = ((PyArrayObject*)raw_input)->dimensions[2];
    }
    else
    {
        n_imgs = 1;
        height = ((PyArrayObject*)raw_input)->dimensions[0];
        width = ((PyArrayObject*)raw_input)->dimensions[1];
    }
    
    raw_input_data  = ((PyArrayObject*)raw_input)->data;
    raw_input_stride = ((PyArrayObject*)raw_input)->strides[((PyArrayObject*)raw_input)->nd - 1];
    
    mask_data = ((PyArrayObject*)mask)->data;
    mask_stride = ((PyArrayObject*)mask)->strides[((PyArrayObject*)mask)->nd - 1];
    
    PyObject * gmf_response = NULL;
    PyObject * gmf_response_angles = NULL;
    if (n_imgs > 1 && par_sigma_scales > 1) 
    {
        npy_intp gmf_response_shape[] = { n_imgs, par_sigma_scales, height, width };
        gmf_response = PyArray_SimpleNew(4, &gmf_response_shape[0], NPY_DOUBLE);
        gmf_response_angles = PyArray_SimpleNew(4, &gmf_response_shape[0], NPY_DOUBLE);
    }
    else if (n_imgs > 1 && par_sigma_scales == 1)
    {
        npy_intp gmf_response_shape[] = { n_imgs, height, width };
        gmf_response = PyArray_SimpleNew(3, &gmf_response_shape[0], NPY_DOUBLE);
        gmf_response_angles = PyArray_SimpleNew(3, &gmf_response_shape[0], NPY_DOUBLE);
    }
    else if (n_imgs == 1 && par_sigma_scales > 1)
    {
        npy_intp gmf_response_shape[] = { par_sigma_scales, height, width };
        gmf_response = PyArray_SimpleNew(3, &gmf_response_shape[0], NPY_DOUBLE);
        gmf_response_angles = PyArray_SimpleNew(3, &gmf_response_shape[0], NPY_DOUBLE);
    }
    else
    {
        npy_intp gmf_response_shape[] = { height, width };
        gmf_response = PyArray_SimpleNew(2, &gmf_response_shape[0], NPY_DOUBLE);
        gmf_response_angles = PyArray_SimpleNew(2, &gmf_response_shape[0], NPY_DOUBLE);
    }
    
    char * gmf_response_data = ((PyArrayObject*)gmf_response)->data;
    char * gmf_response_angles_data = ((PyArrayObject*)gmf_response_angles)->data;
    npy_intp gmf_response_stride = ((PyArrayObject*)gmf_response)->strides[((PyArrayObject*)gmf_response)->nd-1];
    npy_intp gmf_response_angles_stride = ((PyArrayObject*)gmf_response_angles)->strides[((PyArrayObject*)gmf_response_angles)->nd-1];
    
    if (n_imgs > 1)
    {
        DEBMSG("Multiscale Gaussian matched filtering over multiple images\n");
        multiscaleGMFilterWithAngles_multipleinputs((double*)raw_input_data, n_imgs, mask_data, (double*)gmf_response_data, (double*)gmf_response_angles_data, height, width, par_T, par_L, (double*)par_sigma_data, par_sigma_scales, par_K, 0);
    }
    else
    {
        DEBMSG("Multiscale Gaussian matched filtering over a single image\n");
        multiscaleGMFilterWithAngles((double*)raw_input_data, mask_data, (double*)gmf_response_data,(double*)gmf_response_angles_data, height, width, par_T, par_L, (double*)par_sigma_data, par_sigma_scales, par_K, 0);
    }
    
    PyObject *gmf_response_tuple = PyTuple_New(2);
    PyTuple_SetItem(gmf_response_tuple, 0, gmf_response);
    PyTuple_SetItem(gmf_response_tuple, 1, gmf_response_angles);
    return gmf_response_tuple;
}
#endif



#ifdef BUILDING_PYTHON_MODULE
static PyMethodDef gmf_methods[] = {
	{ "gmfFilter", gmfFilter, METH_VARARGS, "applies the Gaussian matched filter to the input image, using the parameters T, L, sigma and K passed, if the parameter sigma is a list, then the multiscale Gaussian matched filter is applied instead." },
	{ "gmfFilterWithAngles", gmfFilterWithAngles, METH_VARARGS, "applies the Gaussian matched filter to the input image, using the parameters T, L, sigma and K passed, if the parameter sigma is a list, then the multiscale Gaussian matched filter is applied instead, the angles of the max response are saved too." },
	{ NULL, NULL, 0, NULL }
};
#endif


#ifdef BUILDING_PYTHON_MODULE
static struct PyModuleDef gmf_moduledef = {
	PyModuleDef_HEAD_INIT,
	"gmf",
	NULL,
	-1,
	gmf_methods,
	NULL,
	NULL,
	NULL,
	NULL
};
#endif


#ifdef BUILDING_PYTHON_MODULE
PyMODINIT_FUNC PyInit_gmf(void)
{
	PyObject *m;
	m = PyModule_Create(&gmf_moduledef);
	if (!m) {
		return NULL;
	}
	import_array();

	return m;
}
#endif
