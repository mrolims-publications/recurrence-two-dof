#ifndef TYPES_H
#define TYPES_H

typedef double (*V_func_t)(const double *, const double *);
typedef double (*T_func_t)(const double *, const double *);
typedef void (*gradV_func_t)(double *, const double *, const double *);
typedef void (*hessV_func_t)(double *, const double *, const double *);
typedef void (*gradT_func_t)(double *, const double *, const double *);
typedef void (*hessT_func_t)(double *, const double *, const double *);

#endif
