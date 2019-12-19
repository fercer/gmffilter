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


double * generateTemplate(double * template_src, const int par_T, const int par_L, unsigned int *kernel_height, unsigned int *kernel_width, const unsigned char untrimmed_kernels)
{	
	unsigned int offset_x, offset_y;

	if (untrimmed_kernels)
	{
		const unsigned int max_dim = (unsigned int)ceil(sqrt(2.0) * ((par_T > par_L) ? par_T : par_L));
		*kernel_height = max_dim;
		*kernel_width = max_dim;

		offset_x = (unsigned int)floor((max_dim - par_T)/2.0);
		offset_y = (unsigned int)floor((max_dim - par_L)/2.0);
	}
	else
	{
		*kernel_height = par_L + 6;
		*kernel_width = par_T + 3;

		offset_x = 1;
		offset_y = 3;
	}

	/* Generate the GMF kernels */
	double sum_GMF_kernel_base = 0.0;
	double * template_src_it = template_src;

	for (unsigned int i = 0; i < par_L; i++)
	{
		for (unsigned int j = 0; j < par_T; j++, template_src_it++)
		{
			sum_GMF_kernel_base += *template_src_it;
		}
	}
	
	const double mean_GMF_kernel_base = sum_GMF_kernel_base / (double)(par_T*par_L);

	double * base_kernel = (double*)calloc(*kernel_height * *kernel_width, sizeof(double));

	double * GMF_kernel_it = base_kernel + offset_y* *kernel_width;
	template_src_it = template_src;
	for (unsigned int i = 0; i < par_L; i++, GMF_kernel_it+=offset_x)
	{
		for (unsigned int j = 0; j < par_T; j++, GMF_kernel_it++, template_src_it++)
		{
			*GMF_kernel_it = (*template_src_it - mean_GMF_kernel_base) / sum_GMF_kernel_base;
		}
	}

	return base_kernel;
}



double * generateGMFTemplate(const int par_T, const int par_L, const double par_sigma, unsigned int *kernel_height, unsigned int *kernel_width, const unsigned char untrimmed_kernels)
{	
	unsigned int offset_x, offset_y;

	if (untrimmed_kernels)
	{
		const unsigned int max_dim = (unsigned int)ceil(sqrt(2.0) * ((par_T > par_L) ? par_T : par_L));
		*kernel_height = max_dim;
		*kernel_width = max_dim;

		offset_x = (unsigned int)floor((max_dim - par_T)/2.0);
		offset_y = (unsigned int)floor((max_dim - par_L)/2.0);

	}
	else
	{
		*kernel_height = par_L + 6;
		*kernel_width = par_T + 3;

		offset_x = 1;
		offset_y = 3;
		
	}

	double * GMF_base_kernel = (double*)calloc(*kernel_height * *kernel_width, sizeof(double));
	
	/* Generate the GMF kernels */
	double curr_x = -floor((double)(par_T) / 2.0);
	double sum_GMF_kernel_base = 0.0;
	double * GMF_kernel_it = GMF_base_kernel + offset_y * *kernel_width + offset_x;

	for (unsigned int j = 0; j < par_T; j++, curr_x+=1.0, GMF_kernel_it++)
	{
		*GMF_kernel_it = 1.0 - exp(-curr_x*curr_x/(2.0*par_sigma*par_sigma));
		sum_GMF_kernel_base += *GMF_kernel_it;
	}
	
	const double mean_GMF_kernel_base = sum_GMF_kernel_base / (double)(par_T);
	GMF_kernel_it = GMF_base_kernel + offset_y * *kernel_width + offset_x;
	for (unsigned int j = 0; j < par_T; j++, GMF_kernel_it++)
	{
		*GMF_kernel_it = (*GMF_kernel_it - mean_GMF_kernel_base) / sum_GMF_kernel_base;
	}

	for (unsigned int i = 1; i < par_L; i++)	
	{
		memcpy(GMF_base_kernel + (i+offset_y)* *kernel_width, GMF_base_kernel + offset_y* *kernel_width, *kernel_width * sizeof(double));
	}

	return GMF_base_kernel;
}



ft_complex ** generateFilterbank(double * base_kernel, const int par_K, const int nearest_2p_dim, const unsigned int kernel_height, const unsigned int kernel_width)
{
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	
	ft_complex ** fft_filter_bank = allocate_ft_complex_pointers(par_K);

	/* Orientate the kernels into K directions */
	/* Transform the zero padded kernels to the frequencies space */
	double *zp_rotated_kernel = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
	double *rotated_kernel = (double*)malloc(kernel_height * kernel_width * sizeof(double));
	
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_filter_bank + k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
		rotateBicubic_impl(base_kernel, rotated_kernel, (double)k * 180.0 / (double)par_K, kernel_height, kernel_width);
		
		for (unsigned int i = 0; i < kernel_height; i++)
		{
			memcpy(zp_rotated_kernel + i * nearest_2p_dim, rotated_kernel + i* kernel_width, kernel_width * sizeof(double));
		}

		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_rotated_kernel, *(fft_filter_bank + k));
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_rotated_kernel, *(fft_filter_bank + k));
		ft_release_forward;
	}

	free(zp_rotated_kernel);
	free(rotated_kernel);
	ft_close;

	return fft_filter_bank;
}



double ** generateFilterbank_space(double * base_kernel, const int par_K, const int nearest_2p_dim, const unsigned int kernel_height, const unsigned int kernel_width)
{	
	double ** filter_bank = (double**)malloc(par_K * sizeof(double*));

	/* Orientate the kernels into K directions */
	/* Transform the zero padded kernels to the frequencies space */
	double *rotated_kernel = (double*)malloc(kernel_height * kernel_width * sizeof(double));
	
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(filter_bank+k) = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
		rotateBicubic_impl(base_kernel, rotated_kernel, (double)k * 180.0 / (double)par_K, kernel_height, kernel_width);
		
		for (unsigned int i = 0; i < kernel_height; i++)
		{
			memcpy(*(filter_bank+k) + i * nearest_2p_dim, rotated_kernel + i* kernel_width, kernel_width * sizeof(double));
		}
	}

	free(rotated_kernel);
	return filter_bank;
}



void applyFilterAngles(ft_complex * fft_img_src, double* img_dst, double* ang_dst, const int nearest_2p_dim, const int height, const int width, const unsigned int s_scales, const int kernel_height, const int kernel_width, const int par_K, ft_complex ** fft_filter_bank)
{	
	/* Offset for the zero padded image */
	const unsigned int offset_y = kernel_height / 2;
	const unsigned int offset_x = kernel_width / 2;

	ft_variables(nearest_2p_dim, nearest_2p_dim);

	/* Arrays for filtering process */
	double * resp_to_GMF_kernels_data = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	double * max_resp = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	double * max_resp_angle = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	ft_complex * fft_convolution_GMF_kernel = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim / 2 + 1));

	/* First Kernel: */
	for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
	{
		const double real_result = ft_real(*(fft_filter_bank)+i) * ft_real(fft_img_src + i) -
			ft_imag(*(fft_filter_bank)+i) * ft_imag(fft_img_src + i);
		const double imag_result = ft_real(*(fft_filter_bank)+i) * ft_imag(fft_img_src + i) +
			ft_real(fft_img_src + i) * ft_imag(*(fft_filter_bank)+i);

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
			const double real_result = ft_real(*(fft_filter_bank + k) + i) * ft_real(fft_img_src + i) -
				ft_imag(*(fft_filter_bank + k) + i) * ft_imag(fft_img_src + i);

			const double imag_result = ft_real(*(fft_filter_bank + k) + i) * ft_imag(fft_img_src + i) + 
				ft_real(fft_img_src + i) * ft_imag(*(fft_filter_bank + k) + i);

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
		for (unsigned int j = 0; j < width; j++)
		{
			*(img_dst + i*width*s_scales + j*s_scales) = *(max_resp + (i + offset_y)*nearest_2p_dim + j + offset_x);
			*(ang_dst + i*width*s_scales + j*s_scales) = *(max_resp_angle + (i + offset_y)*nearest_2p_dim + j + offset_x);
		}
	}
	
	free(max_resp);
	free(max_resp_angle);
}


void multiscaleFilterAngles(double * raw_input, double * output, double * angles, const unsigned int n_inputs, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char untrimmed_kernels, double ** template_src)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	ft_complex** fft_img_src = allocate_ft_complex_pointers(n_inputs);
	for (unsigned int i = 0; i < n_inputs; i++)
	{
		for (unsigned int y = 0; y < height; y++)
		{
			memcpy(zp_img + y*nearest_2p_dim, raw_input + i*height*width + y*width, width * sizeof(double));
		}

		/* Perform the Fourier Transform: */
		*(fft_img_src + i) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));

		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, *(fft_img_src+i));
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, *(fft_img_src+i));
		ft_release_forward;
	}
	
	free(zp_img);
	
	// Compute the filter kernels:
	unsigned int kernel_height, kernel_width;
	ft_complex ** fft_filter_bank;
	double * base_kernel = NULL;
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		if (template_src)
		{
			base_kernel = generateTemplate(template_src + par_T*par_L*s, par_T, par_L, &kernel_height, &kernel_width, untrimmed_kernels);
			fft_filter_bank = generateFilterbank(base_kernel, par_K, nearest_2p_dim, kernel_height, kernel_width);
			free(base_kernel);
		}
		else
		{
			base_kernel = generateGMFTemplate(par_T, par_L, *(par_sigma + s), &kernel_height, &kernel_width, untrimmed_kernels);
			fft_filter_bank = generateFilterbank(base_kernel, par_K, nearest_2p_dim, kernel_height, kernel_width);
			free(base_kernel);
		}

		for (unsigned int i = 0; i < n_inputs; i++)
		{
			// Apply the single-scale filter:
			applyFilterAngles(*(fft_img_src + i), output + i*sigma_scales*width*height + s*height*width, angles + i*sigma_scales*width*height + s*height*width, nearest_2p_dim, height, width, sigma_scales, kernel_height, kernel_width, par_K, fft_filter_bank);
		}

		// Free all memory used for image the filtering:
		for (unsigned int k = 0; k < par_K; k++)
		{
			deallocate_ft_complex(*(fft_filter_bank + k));
		}
		deallocate_ft_complex(fft_filter_bank);
	}
	
	// Free all memory used for image the filtering:
	for (unsigned int i = 0; i < n_inputs; i++)
	{
		deallocate_ft_complex(*(fft_img_src + i));
	}
	deallocate_ft_complex(fft_img_src);
	ft_close;
}



void applyFilter(ft_complex * fft_img_src, double* img_dst, const int nearest_2p_dim, const int height, const int width, const int kernel_height, const int kernel_width, const int par_K, ft_complex ** fft_filter_bank)
{
	/* Offset for the zero padded image */
	const unsigned int offset_y = kernel_height / 2;
	const unsigned int offset_x = kernel_width / 2;

	ft_variables(nearest_2p_dim, nearest_2p_dim);

	/* Arrays for filtering process */
	double * resp_to_GMF_kernels_data = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	double * max_resp = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	ft_complex * fft_convolution_GMF_kernel = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim / 2 + 1));

	/* First Kernel: */
	for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
	{
		const double real_result = ft_real(*(fft_filter_bank)+i) * ft_real(fft_img_src + i) -
			ft_imag(*(fft_filter_bank)+i) * ft_imag(fft_img_src + i);
		const double imag_result = ft_real(*(fft_filter_bank)+i) * ft_imag(fft_img_src + i) +
			ft_real(fft_img_src + i) * ft_imag(*(fft_filter_bank)+i);

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
			const double real_result = ft_real(*(fft_filter_bank + k) + i) * ft_real(fft_img_src + i) -
				ft_imag(*(fft_filter_bank + k) + i) * ft_imag(fft_img_src + i);

			const double imag_result = ft_real(*(fft_filter_bank + k) + i) * ft_imag(fft_img_src + i) + 
				ft_real(fft_img_src + i) * ft_imag(*(fft_filter_bank + k) + i);

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
		for (unsigned int j = 0; j < width; j++)
		{
			*(img_dst + i*width + j) = *(max_resp + (i + offset_y)*nearest_2p_dim + j + offset_x);
		}
	}
	
	free(max_resp);
}


void multiscaleFilter(double * raw_input, double * output, const unsigned int n_inputs, const unsigned int height, const unsigned width, const unsigned int par_T, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char untrimmed_kernels, double ** template_src)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	ft_complex** fft_img_src = allocate_ft_complex_pointers(n_inputs);
	for (unsigned int i = 0; i < n_inputs; i++)
	{
		for (unsigned int y = 0; y < height; y++)
		{
			memcpy(zp_img + y*nearest_2p_dim, raw_input + i*width*height + y*width, width * sizeof(double));
		}

		/* Perform the Fourier Transform: */
		*(fft_img_src + i) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));

		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, *(fft_img_src+i));
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, *(fft_img_src+i));
		ft_release_forward;
	}
	
	free(zp_img);
	
	// Compute the filter kernels:
	unsigned int kernel_height, kernel_width;
	ft_complex ** fft_filter_bank;
	double * base_kernel = NULL;
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		if (template_src)
		{
			base_kernel = generateTemplate(template_src + s*par_T*par_L, par_T, par_L, &kernel_height, &kernel_width, untrimmed_kernels);
			fft_filter_bank = generateFilterbank(base_kernel, par_K, nearest_2p_dim, kernel_height, kernel_width);
			free(base_kernel);
		}
		else
		{
			base_kernel = generateGMFTemplate(par_T, par_L, *(par_sigma + s), &kernel_height, &kernel_width, untrimmed_kernels);
			fft_filter_bank = generateFilterbank(base_kernel, par_K, nearest_2p_dim, kernel_height, kernel_width);
			free(base_kernel);
		}

		for (unsigned int i = 0; i < n_inputs; i++)
		{
			// Apply the single-scale filter:
			applyFilter(*(fft_img_src + i), output + i*sigma_scales*width*height + s*height*width, nearest_2p_dim, height, width, kernel_height, kernel_width, par_K, fft_filter_bank);
		}

		// Free all memory used for image filtering:
		for (unsigned int k = 0; k < par_K; k++)
		{
			deallocate_ft_complex(*(fft_filter_bank + k));
		}
		deallocate_ft_complex(fft_filter_bank);
	}
	
	// Free all memory used for image filtering:
	for (unsigned int i = 0; i < n_inputs; i++)
	{
		deallocate_ft_complex(*(fft_img_src + i));
	}
	deallocate_ft_complex(fft_img_src);
	ft_close;
}





#ifdef BUILDING_PYTHON_MODULE
//static PyObject* gmfFilter(PyObject *self, PyObject *args, PyObject *kwargs)
static PyObject* gmfFilter(PyObject *self, PyObject *args)
{
    PyArrayObject *raw_input;
    char *raw_input_data = NULL;
    npy_intp raw_input_stride, n_imgs = 1;
    npy_intp height, width;
    
    PyArrayObject *multiscale_par_sigma = NULL;
    npy_intp par_sigma_stride, par_sigma = 1;
    
    unsigned int par_T;
    unsigned int par_L;
    unsigned int par_K;
    /*
	static char *kwlist[] = {
        "untrimmed",
        "templates_src",
        NULL
    };
	*/

    //if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!IIO!IO!", kw_list, &PyArray_Type, &raw_input, &par_T, &par_L, &PyArray_Type, &multiscale_par_sigma, &par_K))
    if (!PyArg_ParseTuple(args, "O!IIO!I", &PyArray_Type, &raw_input, &par_T, &par_L, &PyArray_Type, &multiscale_par_sigma, &par_K))
    {
        return NULL;
    }
    
    double * par_sigma_data = (double*)((PyArrayObject*)multiscale_par_sigma)->data;
    par_sigma = ((PyArrayObject*)multiscale_par_sigma)->dimensions[0];
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
    
    raw_input_data  = (double*)((PyArrayObject*)raw_input)->data;
    raw_input_stride = ((PyArrayObject*)raw_input)->strides[((PyArrayObject*)raw_input)->nd - 1];
    
    npy_intp gmf_response_shape[] = { n_imgs, par_sigma, height, width };      
    PyObject * gmf_response = PyArray_SimpleNew(4, &gmf_response_shape[0], NPY_DOUBLE);
    
    double * gmf_response_data = (double*)((PyArrayObject*)gmf_response)->data;
    multiscaleFilter(raw_input_data, gmf_response_data, n_imgs, height, width, par_T, par_L, par_sigma_data, par_sigma, par_K, 1, NULL);
    
    return gmf_response;
}
#endif


#ifdef BUILDING_PYTHON_MODULE
static PyObject* gmfFilterBank(PyObject *self, PyObject *args)
{
    PyArrayObject *raw_input;
    char *raw_input_data = NULL;
    npy_intp raw_input_stride, n_imgs = 1;
    npy_intp height, width;
    
    PyArrayObject *multiscale_par_sigma = NULL;
    npy_intp par_sigma_stride, par_sigma = 1;
    
    unsigned int par_T;
    unsigned int par_L;
    unsigned int par_K;
    
    DEBMSG("Filtering by Gaussian matched filters");
    if (!PyArg_ParseTuple(args, "IIO!I", &par_T, &par_L, &PyArray_Type, &multiscale_par_sigma, &par_K))
    {
        return NULL;
    }
    
    double * par_sigma_data = (double*)((PyArrayObject*)multiscale_par_sigma)->data;
    par_sigma = ((PyArrayObject*)multiscale_par_sigma)->dimensions[0];
    par_sigma_stride = ((PyArrayObject*)multiscale_par_sigma)->strides[((PyArrayObject*)multiscale_par_sigma)->nd-1];
    
	DEBNUMMSG("Sigma: %f", *par_sigma_data);
	DEBNUMMSG(", T: %i", par_T);
	DEBNUMMSG(", L:%i", par_L);
	DEBNUMMSG(", K:%i\n", par_K);

	/* Get the nearest 2-based dimension: */
	const double max_dim = 300.0;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);

	unsigned int kernel_height, kernel_width;
	double * base_kernel = generateGMFTemplate(par_T, par_L, par_sigma, &kernel_height, &kernel_width, 0);
	double **filter_bank = generateFilterbank_space(base_kernel, par_K, nearest_2p_dim, kernel_height, kernel_width);
	free(base_kernel);

    npy_intp gmf_response_shape[] = {par_K, nearest_2p_dim, nearest_2p_dim};
    PyArrayObject * gmf_response = PyArray_SimpleNew(3, &gmf_response_shape[0], NPY_DOUBLE);

	base_kernel = (double*)malloc(nearest_2p_dim*nearest_2p_dim*sizeof(double));

	for (unsigned int k = 0; k < par_K; k++)
	{		
		memcpy((double*)(gmf_response->data) + nearest_2p_dim*nearest_2p_dim*k, *(filter_bank+k), nearest_2p_dim*nearest_2p_dim*sizeof(double));
		free(*(filter_bank + k));
	}
	free(filter_bank);

    return (PyObject*)gmf_response;
}
#endif

#ifdef BUILDING_PYTHON_MODULE
static PyMethodDef gmf_methods[] = {
	{ "gmfFilterBank", gmfFilterBank, METH_VARARGS, "Generate Gaussian matched filter bank." },
	{ "gmfFilter", gmfFilter, METH_VARARGS, "applies the Gaussian matched filter to the input image, using the parameters T, L, sigma and K passed, if the parameter sigma is a list, then the multiscale Gaussian matched filter is applied instead." },
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
