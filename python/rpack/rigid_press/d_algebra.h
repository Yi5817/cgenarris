#ifndef D_ALGEBRA_H
#define D_ALGEBRA_H


void copy_mat3b3_mat3b3(double a[3][3], double b[3][3]);
void inverse_mat3b3(double Tinv[3][3], double T[3][3]);
void mat3b3_transpose(double b_trans[3][3], double b[3][3]);
void vector3_subtract(double a[3], double b[3], double diff[3]);
double det_mat3b3(double a[3][3]);
int hall_number_from_spg(int spg);
void vector3_frac(double a[3]);


static inline void vector3_mat3b3_multiply(double a[3][3], double b[3], double c[3])
{
    double temp[3];
    temp[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
    temp[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
    temp[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
    c[0]  = temp[0];
    c[1]  = temp[1];
    c[2]  = temp[2];
    return;
}

static inline void copy_vector3_vector3(double a[3], double b[3])
{
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
    return;
}

static inline void vector3_intmat3b3_multiply(int a[3][3], double b[3], double c[3])
{
    double temp[3];
    for (int i = 0; i < 3; i++)
        temp[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2];

    c[0] = temp[0];
    c[1] = temp[1];
    c[2] = temp[2];
    return;
}

static inline void vector3_add(double a[3], double b[3], double sum[3])
{
    sum[0] = a[0] + b[0];
    sum[1] = a[1] + b[1];
    sum[2] = a[2] + b[2];
    return;
}

static inline void mat3b3_add(double a[3][3], double b[3][3], double sum[3][3])
{
    sum[0][0] = a[0][0] + b[0][0];
    sum[0][1] = a[0][1] + b[0][1];
    sum[0][2] = a[0][2] + b[0][2];

    sum[1][0] = a[1][0] + b[1][0];
    sum[1][1] = a[1][1] + b[1][1];
    sum[1][2] = a[1][2] + b[1][2];

    sum[2][0] = a[2][0] + b[2][0];
    sum[2][1] = a[2][1] + b[2][1];
    sum[2][2] = a[2][2] + b[2][2];
}

static inline void mat3b3_mat3b3_multiply(double a[3][3], double b[3][3], double c[3][3])
{
    double temp[3][3];
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            temp[i][j] = 0;
            for (int k = 0; k < 3; k++)
                temp[i][j] += a[i][k]*b[k][j];
        }
    }

    copy_mat3b3_mat3b3(c,temp);
    return;
}

#endif
