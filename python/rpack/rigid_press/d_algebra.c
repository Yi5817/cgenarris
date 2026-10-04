#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "d_algebra.h"

#define PI 3.141592653

void rotation_mat_around_axis(double rot[3][3], double axis[3], double psi)
{
    double cosp = cos(psi);
    double vercosp = 1 - cosp; //Thanks S.L. Loney
    double sinp = sin(psi);

    //using Rodrigues rotation formula
    rot[0][0] = cosp + axis[0]*axis[0]*vercosp;
    rot[0][1] = axis[0]*axis[1]*vercosp - axis[2]*sinp;
    rot[0][2] = axis[0]*axis[2]*vercosp + axis[1]*sinp;

    rot[1][0] = axis[0]*axis[1]*vercosp + axis[2]*sinp;
    rot[1][1] = cosp + axis[1]*axis[1]*vercosp;
    rot[1][2] = axis[1]*axis[2]*vercosp - axis[0]*sinp;

    rot[2][0] = axis[0]*axis[2]*vercosp - axis[1]*sinp;
    rot[2][1] = axis[1]*axis[2]*vercosp + axis[0]*sinp;
    rot[2][2] = cosp + axis[2]*axis[2]*vercosp;
}

void print_mat3b3(double mat[3][3])
{
    for(int i = 0; i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
            printf("%f   ", mat[i][j]);
        printf("\n");
    }
}

void print_vec3(double vec[3])
{
    for(int i = 0; i < 3; i++)
        printf("%f   ", vec[i]);
    printf("\n");
}

double det_mat3b3(double a[3][3])
{
     return a[0][0] * ((a[1][1]*a[2][2]) - (a[2][1]*a[1][2])) -a[0][1]
     * (a[1][0] * a[2][2] - a[2][0] * a[1][2]) + a[0][2] * (a[1][0]
     * a[2][1] - a[2][0] * a[1][1]);
}


void inverse_mat3b3(double Tinv[3][3], double T[3][3])
{
    double det = det_mat3b3(T);

    if(det == 0)
       {printf("error: det = 0, matrix has no inverse \n"); print_mat3b3(T);}

    for(int i = 0;i < 3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            Tinv[j][i]=((T[(i+1)%3][(j+1)%3] * T[(i+2)%3][(j+2)%3]) -
                (T[(i+1)%3][(j+2)%3]*T[(i+2)%3][(j+1)%3]))/ det;
        }
    }
}

int hall_number_from_spg(int spg)
{
        static int spg_to_hall_number[230] =
    {
    1,   2,   3,   6,   9,  18,  21,  30,  39,  57,
   60,  63,  72,  81,  90, 108, 109, 112, 115, 116,
  119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
  155, 161, 164, 170, 173, 176, 182, 185, 191, 197,
  203, 209, 212, 215, 218, 221, 227, 228, 230, 233,
  239, 245, 251, 257, 263, 266, 269, 275, 278, 284,
  290, 292, 298, 304, 310, 313, 316, 322, 334, 335,
  337, 338, 341, 343, 349, 350, 351, 352, 353, 354,
  355, 356, 357, 358, 359, 361, 363, 364, 366, 367,
  368, 369, 370, 371, 372, 373, 374, 375, 376, 377,
  378, 379, 380, 381, 382, 383, 384, 385, 386, 387,
  388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
  398, 399, 400, 401, 402, 404, 406, 407, 408, 410,
  412, 413, 414, 416, 418, 419, 420, 422, 424, 425,
  426, 428, 430, 431, 432, 433, 435, 436, 438, 439,
  440, 441, 442, 443, 444, 446, 447, 448, 449, 450,
  452, 454, 455, 456, 457, 458, 460, 462, 463, 464,
  465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
  475, 476, 477, 478, 479, 480, 481, 482, 483, 484,
  485, 486, 487, 488, 489, 490, 491, 492, 493, 494,
  495, 497, 498, 500, 501, 502, 503, 504, 505, 506,
  507, 508, 509, 510, 511, 512, 513, 514, 515, 516,
  517, 518, 520, 521, 523, 524, 525, 527, 529, 530,
    };

    return spg_to_hall_number[spg-1];

}

void copy_mat3b3_mat3b3(double a[3][3], double b[3][3])
{
    a[0][0] = b[0][0];
    a[0][1] = b[0][1];
    a[0][2] = b[0][2];

    a[1][0] = b[1][0];
    a[1][1] = b[1][1];
    a[1][2] = b[1][2];

    a[2][0] = b[2][0];
    a[2][1] = b[2][1];
    a[2][2] = b[2][2];

    return;
}


 void mat3b3_transpose(double b_trans[3][3], double b[3][3])
{
    double temp[3][3];

    temp[0][0] = b[0][0];
    temp[1][1] = b[1][1];
    temp[2][2] = b[2][2];

    temp[1][0] = b[0][1];
    temp[2][0] = b[0][2];
    temp[2][1] = b[1][2];

    temp[0][1] = b[1][0];
    temp[0][2] = b[2][0];
    temp[1][2] = b[2][1];

    copy_mat3b3_mat3b3(b_trans, temp);
}

void copy_vector3_mat3b3(double a[3], double b[3][3],int index)
{
    a[0] = b[index][0];
    a[1] = b[index][1];
    a[2] = b[index][2];

}

void vector3_inverse(double a[3])
{
    a[0] = -a[0];
    a[1] = -a[1];
    a[2] = -a[2];
}

double vector3_norm(double a[3])
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

void vector3_subtract(double a[3], double b[3], double diff[3])
{
    diff[0] = a[0] - b[0];
    diff[1] = a[1] - b[1];
    diff[2] = a[2] - b[2];
    return;
}

void cross_vector3_vector3(double cross[3], double a[3], double b[3])
{
    cross[0] = a[1]*b[2] - a[2]*b[1];
    cross[1] = a[2]*b[0] - a[0]*b[2];
    cross[2] = a[0]*b[1] - a[1]*b[0];

    return;
}

void normalise_vector3(double a[3])
{
    double norm;
    norm = vector3_norm(a);
    a[0] /= norm;
    a[1] /= norm;
    a[2] /= norm;
}

double dot_vector3_vector3(double a[3], double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

int check_vec3_isNull(double a[3], double tol)
{
    if (fabs(a[0]) < tol &&
        fabs(a[1]) < tol &&
        fabs(a[2]) < tol  )
        return 1;
    else
        return 0;
}

double cart_dist(double p1[3], double p2[3])
{
        double c[3];
        vector3_subtract(p1,p2,c);
        return vector3_norm(c);
}


int are_equal_doubles(double a, double b, double ftol)
{
    if(fabs(a - b) < ftol)
        return 1;
    return 0;
}

/*
Get fractional part of the vector.
*/
void vector3_frac(double a[3])
{
    int b[3];
    double frac;
    for(int i = 0; i < 3; i++)
    {
	b[i] = (int) a[i];

	if(a[i] >= 0)
	{
	    frac = a[i] - b[i];
	}
	else
	{
	    frac = a[i] - b[i] + 1;
	}

	a[i] = frac;
    }
}


// get integer part of a vector
void vector3_int(double a[3])
{
    int b[3];
    for(int i = 0; i < 3; i++)
    {
        b[i] = (int)a[i];

        if(a[i] < 0)
            b[i] -= 1;

        a[i] = b[i];
    }
}

