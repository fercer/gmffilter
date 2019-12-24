#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmf.h>

void savePGM(const char * filename, const int height, const int width, const double * img_src);

int main(int argc, char * argv[])
{    
    if (argc < 2)
    {
        printf("\nExample: Filter an input image using the GMF filter\n");
        printf("First argument is the path to a .PGM image in raw format\n");
        printf("Second to fourth arguments are the parameters for the GMF filter\n");
        return -1;
    }
 
    int par_K = 12;
    int par_L = 9;
    int par_T = 13;
    double par_sigma = 2.0;

    FILE * fp_img = fopen(argv[1], "r");
    if (!fp_img)
    {
        printf("<<Error: The file \"%s\" could not be loaded>>\n", argv[1]);
        return -1;
    }
    
    char magic_number[3];
    fgets(magic_number, 3, fp_img);
    if (strcmp(magic_number, "P5") && strcmp(magic_number, "P2"))
    {
        printf("<<Error: The input file is not of the required format, its format is: \"%s\" instead of \"P5/P2\">>\n", magic_number);
        fclose(fp_img);
        return -1;
    }
    fgetc(fp_img);
    
    char commentary[64];
    fgets(commentary, 64, fp_img);
    printf("Commentary: %s\n", commentary);
    int height, width;
    if (*commentary == '#')
    {
        fscanf(fp_img, "%i", &height);
        fscanf(fp_img, "%i", &width);
    }
    else
    {
        height = atoi(commentary);
        char * width_ptr = strchr(commentary, ' ') + 1;
        width = atoi(width_ptr);
    }
    
    int max_value;
    fscanf(fp_img,"%i", &max_value);
    
    printf("The image has dimension: %ix%i pixels, with maximum value of: %i\n", height, width, max_value);
        
    int read_value;
    char * mask = (char*) malloc(height * width * sizeof(char));
    double * img = (double*) malloc(height * width * sizeof(double));
    double * resp = (double*) malloc(height * width * sizeof(double));
    double * ang_resp = (double*) malloc(height * width * sizeof(double));

    if (magic_number[1] == '2')
    {
        printf("Reading the image from an ascii file\n");
        for (unsigned int xy = 0; xy < height*width; xy++)
        {
            fscanf(fp_img, "%i", &read_value);
            *(img+xy) = (double)read_value / 255.0;
            *(mask+xy) = (char)1;
        }
    }
    else
    {
        printf("Reading the image from a raw file\n");
        for (unsigned int xy = 0; xy < height*width; xy++)
        {
            read_value = fgetc(fp_img);
            *(img+xy) = (double)read_value / 255.0;
            *(mask+xy) = (char)1;
        }
    }
    
    fclose(fp_img);

    savePGM("loaded_image.pgm", height, width, img);
    
    printf("Image loaded (first, middle and last values are: %.16f, %.16f, %.16f)...\n", *img, *(img+(150-1)*300+150-1), *(img + height*width-1));
    
    singleScaleGMFilterWithAngles(img, mask, resp, ang_resp, height, width, par_T, par_L, par_sigma, par_K, 0);

    free(img);
    free(mask);
    
    fp_img = fopen("gmf_resp.bin", "wb");
    
    fwrite(&height, sizeof(int), 1, fp_img);
    fwrite(&width, sizeof(int), 1, fp_img);
    fwrite(resp, sizeof(double), height*width, fp_img);
    fwrite(ang_resp, sizeof(double), height*width, fp_img);
    
    fclose(fp_img);

    printf("Dumping the image to an ascii file\n");

    savePGM("resp.pgm", height, width, resp);
    savePGM("ang_resp.pgm", height, width, ang_resp);

    free(resp);   
    free(ang_resp);   
    
    printf("Example run successfully ...\n");
   
    return 0;
}






void savePGM(const char * filename, const int height, const int width, const double * img_src)
{
    double max_value = -MY_INF;
    double min_value = MY_INF;
    
    for (unsigned int xy = 0; xy < height*width; xy++)
    {
        if (max_value < *(img_src + xy))
        {
            max_value = *(img_src + xy);
        }
        
        if (min_value > *(img_src + xy))
        {
            min_value = *(img_src + xy);
        }
    }
    
    FILE * fp = fopen(filename, "w");
    fprintf(fp, "P2\n# Created by FerCer\n%i %i\n255\n", height, width);
    
    for (unsigned int x = 0; x < width; x++)
    {
        for (unsigned int y = 0; y < height; y++)
        {
            fprintf(fp, "%i\n", (int)(255.0 * (*(img_src + x + y*width) - min_value) / (max_value - min_value)));
        }
    }
    
    fclose(fp);
}