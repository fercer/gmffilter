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


void generateTemplate(double * template_filter, const unsigned int par_L, const double par_sigma, const double par_theta, unsigned int nearest_2p_dim)
{
	/* Generate the DOG kernels */
	double y, x, u, v;
	double dog_template, step_upper, step_lower;
	const double ctheta = cos(par_theta), stheta = sin(par_theta);
	/* T in GMF is 6 * sigma, and the area (volume) under the filter curve is T * L */
	const double alpha = 1.0 + 1.0/(6.0*par_sigma*par_L);
	y = -(double)nearest_2p_dim / 2.0;

	for (unsigned int i = 0; i < nearest_2p_dim; i++, y+=1.0)
	{
		x = -(double)nearest_2p_dim / 2.0;
		for (unsigned int j = 0; j < nearest_2p_dim; j++, x+=1.0)
		{
			u = x * ctheta - y * stheta;
			v = y * ctheta + x * stheta;

			step_upper = 1.0/(1.0 + exp(- v - (double)par_L/2.0));
			step_lower = 1.0/(1.0 + exp(- v + (double)par_L/2.0));
			dog_template = exp(-0.5*u*u/(par_sigma * par_sigma)) - alpha * exp(-0.5*alpha*alpha*u*u/(par_sigma*par_sigma));

			*(template_filter + i * nearest_2p_dim + j) = (step_upper - step_lower) * dog_template;
		}
	}
}



ft_complex ** generateFilterbank(const unsigned int par_L, const double par_sigma, const int par_K, const int nearest_2p_dim)
{
	ft_variables(nearest_2p_dim, nearest_2p_dim);
	
	ft_complex ** fft_filter_bank = allocate_ft_complex_pointers(par_K);

	/* Orientate the kernels into K directions */
	/* Transform the zero padded kernels to the frequencies space */
	double *template_filter = (double*)calloc(nearest_2p_dim * nearest_2p_dim, sizeof(double));
	
	for (unsigned int k = 0; k < par_K; k++)
	{
		*(fft_filter_bank + k) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));
		generateTemplate(template_filter, par_L, par_sigma, (double)k/180.0 * MY_PI, nearest_2p_dim);

		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, template_filter, *(fft_filter_bank + k));
		ft_forward(nearest_2p_dim, nearest_2p_dim, template_filter, *(fft_filter_bank + k));
		ft_release_forward;
	}

	free(template_filter);
	ft_close;

	return fft_filter_bank;
}


void applyFilterAngles(ft_complex * fft_img_src, double* img_dst, double* ang_dst, const int nearest_2p_dim, const int height, const int width, const unsigned int s_scales, const int par_K, ft_complex ** fft_filter_bank)
{	
	/* Offset for the zero padded image */
	const unsigned int offset_y = (nearest_2p_dim - height) / 2;
	const unsigned int offset_x = (nearest_2p_dim - width) / 2;

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


void multiscaleFilterAngles(double * raw_input, double * output, double * angles, const unsigned int n_imgs, const unsigned int height, const unsigned width, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	ft_complex** fft_img_src = allocate_ft_complex_pointers(n_imgs);
	for (unsigned int i = 0; i < n_imgs; i++)
	{
		for (unsigned int y = 0; y < height/2; y++)
		{
			memcpy(zp_img + (y+nearest_2p_dim/2 + height/2)*nearest_2p_dim + nearest_2p_dim/2 + width/2, raw_input + i*height*width + y*width, width / 2 * sizeof(double));
			memcpy(zp_img + (y+nearest_2p_dim/2 + height/2)*nearest_2p_dim, raw_input + i*height*width + y*width + width/2, width / 2 * sizeof(double));
		}

		for (unsigned int y = height/2; y < height; y++)
		{
			memcpy(zp_img + y*nearest_2p_dim + nearest_2p_dim/2 + width/2, raw_input + i*height*width + y*width, width / 2 * sizeof(double));
			memcpy(zp_img + y*nearest_2p_dim, raw_input + i*height*width + y*width + width/2, width / 2 * sizeof(double));
		}

		/* Perform the Fourier Transform: */
		*(fft_img_src + i) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));

		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, *(fft_img_src+i));
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, *(fft_img_src+i));
		ft_release_forward;
	}
	
	free(zp_img);
	
	// Compute the filter kernels:
	ft_complex ** fft_filter_bank;
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		fft_filter_bank = generateFilterbank(par_L, *(par_sigma + s), par_K, nearest_2p_dim);

		for (unsigned int i = 0; i < n_imgs; i++)
		{
			// Apply the single-scale filter:
			applyFilterAngles(*(fft_img_src + i), output + i*sigma_scales*width*height + s*height*width,  angles + i*sigma_scales*width*height + s*height*width, nearest_2p_dim, height, width, sigma_scales, par_K, fft_filter_bank);
		}

		// Free all memory used for image the filtering:
		for (unsigned int k = 0; k < par_K; k++)
		{
			deallocate_ft_complex(*(fft_filter_bank + k));
		}
		deallocate_ft_complex(fft_filter_bank);
	}
	
	// Free all memory used for image the filtering:
	for (unsigned int i = 0; i < n_imgs; i++)
	{
		deallocate_ft_complex(*(fft_img_src + i));
	}
	deallocate_ft_complex(fft_img_src);
	ft_close;
}



void applyFilter(ft_complex * fft_img_src, double* img_dst, const int nearest_2p_dim, const int height, const int width, const int par_K, ft_complex ** fft_filter_bank)
{
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
			*(img_dst + i*width + j) = *(max_resp + (i + (nearest_2p_dim-height)/2)*nearest_2p_dim + j + (nearest_2p_dim-width)/2);
		}
	}
	
	free(max_resp);
}



void applyTemplates(ft_complex * fft_img_src, double* img_dst, const int nearest_2p_dim, const int height, const int width, const int par_K, ft_complex ** fft_filter_bank)
{
	/* Offset for the zero padded image */
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	/* Arrays for filtering process */
	double * resp_to_GMF_kernels_data = (double*)malloc(nearest_2p_dim * nearest_2p_dim * sizeof(double));
	ft_complex * fft_convolution_GMF_kernel = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim / 2 + 1));
	
	for (unsigned int k = 0; k < par_K; k++)
	{
		for (unsigned int i = 0; i < nearest_2p_dim*(nearest_2p_dim / 2 + 1); i++)
		{
			const double real_result = ft_real(*(fft_filter_bank + k) + i) * ft_real(fft_img_src + i) -
				ft_imag(*(fft_filter_bank + k) + i) * ft_imag(fft_img_src + i);

			const double imag_result = ft_real(*(fft_filter_bank + k) + i) * ft_imag(fft_img_src + i) + 
				ft_real(fft_img_src + i) * ft_imag(*(fft_filter_bank + k) + i);

			ft_real_assign(fft_convolution_GMF_kernel + i) = real_result;
			ft_imag_assign(fft_convolution_GMF_kernel + i) = imag_result;
		}

		ft_backward_setup(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
		ft_backward(nearest_2p_dim, nearest_2p_dim, fft_convolution_GMF_kernel, resp_to_GMF_kernels_data);
		ft_release_backward;

		for (unsigned int i = 0; i < height; i++)
		{
			for (unsigned int j = 0; j < width; j++)
			{
				*(img_dst + k*height*width + i*width + j) = *(resp_to_GMF_kernels_data + (i + (nearest_2p_dim-height)/2)*nearest_2p_dim + j + (nearest_2p_dim-width)/2)  / (double)(nearest_2p_dim*nearest_2p_dim);
			}
		}
	}
	deallocate_ft_complex(fft_convolution_GMF_kernel);
	
	free(resp_to_GMF_kernels_data);
	ft_close;
}


void multiscaleFilter(double * raw_input, double * output, const unsigned int n_imgs, const unsigned int height, const unsigned width, const unsigned int par_L, double * par_sigma, const unsigned int sigma_scales, const unsigned int par_K, const unsigned char compute_max)
{
	/* Get the nearest 2-based dimension: */
	const double max_dim = (height > width) ? (double)height : (double)width;
	const double nearest_power = floor(log2(max_dim)) + 1.0;
	const unsigned int nearest_2p_dim = (unsigned int)pow(2.0, nearest_power);
	
	/* Zero pad the raw input image: */
	double * zp_img = (double*)calloc(nearest_2p_dim*nearest_2p_dim, sizeof(double));
	ft_variables(nearest_2p_dim, nearest_2p_dim);

	ft_complex** fft_img_src = allocate_ft_complex_pointers(n_imgs);
	for (unsigned int i = 0; i < n_imgs; i++)
	{
		for (unsigned int y = 0; y < height/2; y++)
		{
			memcpy(zp_img + (y + nearest_2p_dim - height/2)*nearest_2p_dim + nearest_2p_dim - width/2, raw_input + i*height*width + y*width, width/2 * sizeof(double));
			memcpy(zp_img + (y + nearest_2p_dim - height/2)*nearest_2p_dim, raw_input + i*height*width + y*width + width/2, width/2 * sizeof(double));
		}

		for (unsigned int y = height/2; y < height; y++)
		{
			memcpy(zp_img + (y - height/2)*nearest_2p_dim + nearest_2p_dim - width/2, raw_input + i*height*width + y*width, width/2 * sizeof(double));
			memcpy(zp_img + (y - height/2)*nearest_2p_dim, raw_input + i*height*width + y*width + width/2, width/2 * sizeof(double));
		}

		FILE * fp_img = fopen("test.bin","wb");
		fwrite(zp_img, sizeof(double), nearest_2p_dim*nearest_2p_dim,fp_img);
		fclose(fp_img);

		/* Perform the Fourier Transform: */
		*(fft_img_src + i) = allocate_ft_complex(nearest_2p_dim * (nearest_2p_dim/2 + 1));

		ft_forward_setup(nearest_2p_dim, nearest_2p_dim, zp_img, *(fft_img_src+i));
		ft_forward(nearest_2p_dim, nearest_2p_dim, zp_img, *(fft_img_src+i));
		ft_release_forward;
	}
	
	free(zp_img);
	
	// Compute the filter kernels:
	ft_complex ** fft_filter_bank;
	for (unsigned int s = 0; s < sigma_scales; s++)
	{
		fft_filter_bank = generateFilterbank(par_L, *(par_sigma + s), par_K, nearest_2p_dim);
		if (compute_max)
		{
			for (unsigned int i = 0; i < n_imgs; i++)
			{
				// Apply the single-scale filter to get the max response:
				applyFilter(*(fft_img_src + i), output + i*sigma_scales*width*height + s*height*width, nearest_2p_dim, height, width, par_K, fft_filter_bank);
			}
		}
		else
		{
			for (unsigned int i = 0; i < n_imgs; i++)
			{
				// Apply the single-scale filter:
				applyTemplates(*(fft_img_src + i), output + i*sigma_scales*par_K*width*height + s*par_K*height*width, nearest_2p_dim, height, width, par_K, fft_filter_bank);
			}
		}

		// Free all memory used for image filtering:
		for (unsigned int k = 0; k < par_K; k++)
		{
			deallocate_ft_complex(*(fft_filter_bank + k));
		}
		deallocate_ft_complex(fft_filter_bank);
	}
	
	// Free all memory used for image filtering:
	for (unsigned int i = 0; i < n_imgs; i++)
	{
		deallocate_ft_complex(*(fft_img_src + i));
	}
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
    npy_intp par_sigma_stride, sigma_scales = 1;
    
    unsigned int par_L;
    unsigned int par_K;
	unsigned char compute_max = 0;

	/** untrimmed kernel indicator and template src is an optional argument, wich is specified after '|' */
    if (!PyArg_ParseTuple(args, "O!IO!I|b", &PyArray_Type, &raw_input, &par_L, &PyArray_Type, &multiscale_par_sigma, &par_K, &compute_max))
    {
        return NULL;
    }
    
    double * par_sigma_data = (double*)((PyArrayObject*)multiscale_par_sigma)->data;
    sigma_scales = ((PyArrayObject*)multiscale_par_sigma)->dimensions[0];
    par_sigma_stride = ((PyArrayObject*)multiscale_par_sigma)->strides[((PyArrayObject*)multiscale_par_sigma)->nd-1];
    
    if (((PyArrayObject*)raw_input)->nd > 3)
    {
        n_imgs = ((PyArrayObject*)raw_input)->dimensions[0];
		/* Do not consider channels */
        height = ((PyArrayObject*)raw_input)->dimensions[2];
        width = ((PyArrayObject*)raw_input)->dimensions[3];
    }        
    else if (((PyArrayObject*)raw_input)->nd > 2)
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
    
    npy_intp gmf_response_shape[] = { n_imgs, sigma_scales * (compute_max ? 1 : par_K), height, width };
    PyObject * gmf_response = PyArray_SimpleNew(4, &gmf_response_shape[0], NPY_DOUBLE);
    
    double * gmf_response_data = (double*)((PyArrayObject*)gmf_response)->data;
	printf("[Filtering image of shape: {%i, %i, %i} with maximum: %i --> %f]\n", n_imgs, height, width, sigma_scales * (compute_max ? 1 : par_K), *gmf_response_data);
	multiscaleFilter(raw_input_data, gmf_response_data, n_imgs, height, width, par_L, par_sigma_data, sigma_scales, par_K, compute_max);
	
    return gmf_response;
}
#endif


#ifdef BUILDING_PYTHON_MODULE
static PyMethodDef gmf_methods[] = {
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
