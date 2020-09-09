#include "helpers.h"
#include "math.h"
#include "stdio.h"

//RGB value counting function to make code for "Blur" more concise
void rgbsum(int height, int width, int i, int j, float *numerator, int x, float sum_rgb[3], RGBTRIPLE image[i][j]);

//RGB value counting function to make code for "Edge" more concise
void Gx_sum(int height, int width, int i, int j, int x, float sum_Gx[3], RGBTRIPLE image[i][j], int G[3][3], int row);
void Gy_sum(int height, int width, int i, int j, int x, float sum_Gy[3], RGBTRIPLE image[i][j], int G[3][3], int row);

// Convert image to grayscale
void grayscale(int height, int width, RGBTRIPLE image[height][width])
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            //Calculate average value of the r,g, and b component of each pixel
            float grey_value = (float)(image[i][j].rgbtBlue + image[i][j].rgbtGreen + image[i][j].rgbtRed) / 3;

            //set each pixel equal to average value to creat a "grey" pixel of appropriate intensity
            image[i][j].rgbtBlue = round(grey_value);
            image[i][j].rgbtGreen = round(grey_value);
            image[i][j].rgbtRed = round(grey_value);
        }
    }
    return;
}

// Reflect image horizontally
void reflect(int height, int width, RGBTRIPLE image[height][width])
{
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < (width / 2); j++)
        {
            RGBTRIPLE temp = image[i][j];
            image[i][j] = image[i][(width - 1) - j];
            image[i][(width - 1) - j] = temp;
        }
    }
    return;
}

// Blur image
void blur(int height, int width, RGBTRIPLE image[height][width])
{
    //RGBTRIPLE sum;
    float sum_rgb[3];
    float numerator;
    RGBTRIPLE holder[height][width];

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            numerator = 0;
            sum_rgb[0] = 0;
            sum_rgb[1] = 0;
            sum_rgb[2] = 0;

            //average top row
            if (i != 0)
            {
                rgbsum(height, width, i, j, &numerator, -1, sum_rgb, image);
            }

            //average bottom row
            if (i != (height - 1))
            {
                rgbsum(height, width, i, j, &numerator, 1, sum_rgb, image);
            }

            //average middle row
            rgbsum(height, width, i, j, &numerator, 0, sum_rgb, image);

            //fill components into RGBTRIPLE holder
            holder[i][j].rgbtBlue = round(sum_rgb[0] / numerator);
            holder[i][j].rgbtGreen = round(sum_rgb[1] / numerator);
            holder[i][j].rgbtRed = round(sum_rgb[2] / numerator);
        }
    }

    //Transfer holder into image
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            image[i][j] = holder[i][j];
        }
    }

    return;
}

void rgbsum(int height, int width, int i, int j, float *numerator, int x, float sum_rgb[], RGBTRIPLE image[height][width])
{
    if (j == 0)
    {
        sum_rgb[0] += (image[i + x][j].rgbtBlue + image[i + x][j + 1].rgbtBlue);
        sum_rgb[1] += (image[i + x][j].rgbtGreen + image[i + x][j + 1].rgbtGreen);
        sum_rgb[2] += (image[i + x][j].rgbtRed + image[i + x][j + 1].rgbtRed);
        *numerator += 2;
    }
    else if (j == (width - 1))
    {
        sum_rgb[0] += (image[i + x][j - 1].rgbtBlue + image[i + x][j].rgbtBlue);
        sum_rgb[1] += (image[i + x][j - 1].rgbtGreen + image[i + x][j].rgbtGreen);
        sum_rgb[2] += (image[i + x][j - 1].rgbtRed + image[i + x][j].rgbtRed);
        *numerator += 2;
    }
    else
    {
        sum_rgb[0] += (image[i + x][j - 1].rgbtBlue + image[i + x][j].rgbtBlue + image[i + x][j + 1].rgbtBlue);
        sum_rgb[1] += (image[i + x][j - 1].rgbtGreen + image[i + x][j].rgbtGreen + image[i + x][j + 1].rgbtGreen);
        sum_rgb[2] += (image[i + x][j - 1].rgbtRed + image[i + x][j].rgbtRed + image[i + x][j + 1].rgbtRed);
        *numerator += 3;
    }
}


// Detect edges
void edges(int height, int width, RGBTRIPLE image[height][width])
{
    //Define Kernels
    int Gx[3][3] = {{-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1}};
    int Gy[3][3] = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
    int G_val[] = {0, 0, 0};
    int row = 0;

    //RGBTRIPLE sum;
    float sum_Gx[3];
    float sum_Gy[3];
    RGBTRIPLE holder[height][width];

    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            sum_Gx[0] = 0, sum_Gx[1] = 0, sum_Gx[2] = 0;
            sum_Gy[0] = 0, sum_Gy[1] = 0, sum_Gy[2] = 0;;

            //Calculate top row sum
            if (i != 0)
            {
                row = 0;
                Gx_sum(height, width, i, j, -1, sum_Gx, image, Gx, row);
                Gy_sum(height, width, i, j, -1, sum_Gy, image, Gy, row);
            }

            //Calculate bottom row sum
            if (i != (height - 1))
            {
                row = 2;
                Gx_sum(height, width, i, j, 1, sum_Gx, image, Gx, row);
                Gy_sum(height, width, i, j, 1, sum_Gy, image, Gy, row);
            }

            //Calculate middle row sum
            row = 1;
            Gx_sum(height, width, i, j, 0, sum_Gx, image, Gx, row);
            Gy_sum(height, width, i, j, 0, sum_Gy, image, Gy, row);

            //Combine Gx and Gy to a single value -> must be a rounded integer capped at 255
            for (int k = 0; k < 3; k++)
            {
                sum_Gx[k] = round(sum_Gx[k]);
                sum_Gy[k] = round(sum_Gy[k]);

                G_val[k] = sum_Gx[k] * sum_Gx[k] + sum_Gy[k] * sum_Gy[k];

                G_val[k] = round(sqrt(G_val[k]));
                //G_val[k] = round(G_val[k]);

                if (G_val[k] > 255)
                {
                    G_val[k] = 255;
                }

            }

            //fill components into RGBTRIPLE holder
            holder[i][j].rgbtBlue = round(G_val[0]);
            holder[i][j].rgbtGreen = round(G_val[1]);
            holder[i][j].rgbtRed = round(G_val[2]);
        }
    }

    //Transfer holder into image
    for (int i = 0; i < height; i++)
    {
        for (int j = 0; j < width; j++)
        {
            image[i][j] = holder[i][j];
        }
    }
    return;
}


void Gx_sum(int height, int width, int i, int j, int a, float sum_Gx[], RGBTRIPLE image[height][width], int G[3][3], int row)
{
    if (j == 0)
    {
        sum_Gx[0] += (image[i + a][j].rgbtBlue * (G[row][1])) + (image[i + a][j + 1].rgbtBlue * (G[row][2]));
        sum_Gx[1] += (image[i + a][j].rgbtGreen * (G[row][1])) + (image[i + a][j + 1].rgbtGreen * (G[row][2]));
        sum_Gx[2] += (image[i + a][j].rgbtRed * (G[row][1])) + (image[i + a][j + 1].rgbtRed * (G[row][2]));
    }
    else if (j == (width - 1))
    {
        sum_Gx[0] += (image[i + a][j - 1].rgbtBlue * (G[row][0])) + (image[i + a][j].rgbtBlue * (G[row][1]));
        sum_Gx[1] += (image[i + a][j - 1].rgbtGreen * (G[row][0])) + (image[i + a][j].rgbtGreen * (G[row][1]));
        sum_Gx[2] += (image[i + a][j - 1].rgbtRed * (G[row][0])) + (image[i + a][j].rgbtRed * (G[row][1]));
    }
    else
    {
        sum_Gx[0] += (image[i + a][j - 1].rgbtBlue * (G[row][0])) + (image[i + a][j].rgbtBlue * (G[row][1])) +
                     (image[i + a][j + 1].rgbtBlue * (G[row][2]));
        sum_Gx[1] += (image[i + a][j - 1].rgbtGreen * (G[row][0])) + (image[i + a][j].rgbtGreen * (G[row][1])) +
                     (image[i + a][j + 1].rgbtGreen * (G[row][2]));
        sum_Gx[2] += (image[i + a][j - 1].rgbtRed * (G[row][0])) + (image[i + a][j].rgbtRed * (G[row][1])) +
                     (image[i + a][j + 1].rgbtRed * (G[row][2]));
    }
}

void Gy_sum(int height, int width, int i, int j, int a, float sum_Gy[], RGBTRIPLE image[height][width], int G[3][3], int row)
{

    if (j == 0)
    {
        sum_Gy[0] += (image[i + a][j].rgbtBlue * (G[row][1])) + (image[i + a][j + 1].rgbtBlue * (G[row][2]));
        sum_Gy[1] += (image[i + a][j].rgbtGreen * (G[row][1])) + (image[i + a][j + 1].rgbtGreen * (G[row][2]));
        sum_Gy[2] += (image[i + a][j].rgbtRed * (G[row][1])) + (image[i + a][j + 1].rgbtRed * (G[row][2]));
    }
    else if (j == (width - 1))
    {
        sum_Gy[0] += (image[i + a][j - 1].rgbtBlue * (G[row][0])) + (image[i + a][j].rgbtBlue * (G[row][1]));
        sum_Gy[1] += (image[i + a][j - 1].rgbtGreen * (G[row][0])) + (image[i + a][j].rgbtGreen * (G[row][1]));
        sum_Gy[2] += (image[i + a][j - 1].rgbtRed * (G[row][0])) + (image[i + a][j].rgbtRed * (G[row][1]));
    }
    else
    {
        sum_Gy[0] += (image[i + a][j - 1].rgbtBlue * (G[row][0])) + (image[i + a][j].rgbtBlue * (G[row][1])) +
                     (image[i + a][j + 1].rgbtBlue * (G[row][2]));
        sum_Gy[1] += (image[i + a][j - 1].rgbtGreen * (G[row][0])) + (image[i + a][j].rgbtGreen * (G[row][1])) +
                     (image[i + a][j + 1].rgbtGreen * (G[row][2]));
        sum_Gy[2] += (image[i + a][j - 1].rgbtRed * (G[row][0])) + (image[i + a][j].rgbtRed * (G[row][1])) +
                     (image[i + a][j + 1].rgbtRed * (G[row][2]));
    }
}