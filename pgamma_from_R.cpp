/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */


#include <math.h>
#include <float.h>
#include <stdio.h>
#include <iostream>


/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
 *  Copyright (C) 2005-10 The R Foundation
 *  Copyright (C) 2006-2015 The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *    #include <Rmath.h>
 *
 *    double pgamma (double x, double alph, double scale,
 *               int lower_tail, int log_p)
 *
 *    double log1pmx    (double x)
 *    double lgamma1p (double a)
 *
 *    double logspace_add (double logx, double logy)
 *    double logspace_sub (double logx, double logy)
 *    double logspace_sum (double* logx, int n)
 *
 *
 *  DESCRIPTION
 *
 *    This function computes the distribution function for the
 *    gamma distribution with shape parameter alph and scale parameter
 *    scale.    This is also known as the incomplete gamma function.
 *    See Abramowitz and Stegun (6.5.1) for example.
 *
 *  NOTES
 *
 *    Complete redesign by Morten Welinder, originally for Gnumeric.
 *    Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
 *
 *  REFERENCES
 *
 */

/*----------- DEBUGGING -------------
 * make CFLAGS='-DDEBUG_p -g'
 * (cd `R-devel RHOME`/src/nmath; gcc -I. -I../../src/include -I../../../R/src/include  -DHAVE_CONFIG_H -fopenmp -DDEBUG_p -g -c ../../../R/src/nmath/pgamma.c -o pgamma.o)
 */

/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define M_LN_SQRT_2PI   0.918938533204672741780329736406    /* log(sqrt(2*pi)) */
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x)))
#define give_log log_p
#define R_D__0    (log_p ? -INFINITY : 0.)
#define R_D__1    (log_p ? 0. : 1.)            /* 1 */
#define R_D_exp(x)    (log_p    ?  (x)     : exp(x))    /* exp(x) */
#define SQR(x) ((x)*(x))
#define R_D_fexp(f,x) (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

#ifdef HAVE_LONG_DOUBLE
#  define LDOUBLE long double
#else
#  define LDOUBLE double
#endif

#ifdef HAVE_NEARYINT
# define R_forceint(x)   nearbyint()
#else
# define R_forceint(x)   round(x)
#endif

#define R_P_bounds_01(x, x_min, x_max)    \
    if(x <= x_min) return R_DT_0;        \

/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]     =~=  -x */
static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxilary in log1pmx() and lgamma1p()
 */

#ifdef HAVE_VISIBILITY_ATTRIBUTE
# define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
# define attribute_hidden
#endif
#define M_LN_SQRT_PId2    0.225791352644727432363097614947
#define     M_2PI   6.28318530717958647692528676655900576
#define M_1_SQRT_2PI    0.398942280401432677939946059934    /* 1/sqrt(2pi) */
#define M_SQRT_32    5.656854249492380195206754896838
#define R_DT_0    (lower_tail ? R_D__0 : R_D__1)        /* 0 */
#define R_DT_1    (lower_tail ? R_D__1 : R_D__0)        /* 1 */

#ifdef HAVE_SINPI
#elif defined HAVE___SINPI
double sinpi(double x) {
    return __sinpi(x);
}
#else
// sin(pi * x)  -- exact when x = k/2  for all integer k
double sinpi(double x) {
#ifdef IEEE_754
    if (ISNAN(x)) return x;
#endif
    if(isinf(x)) return std::numeric_limits<double>::quiet_NaN();

    x = fmod(x, 2.); // sin(pi(x + 2k)) == sin(pi x)  for all integer k
    // map (-2,2) --> (-1,1] :
    if(x <= -1) x += 2.; else if (x > 1.) x -= 2.;
    if(x == 0. || x == 1.) return 0.;
    if(x ==  0.5)    return  1.;
    if(x == -0.5)    return -1.;
    // otherwise
    return sin(M_PI * x);
}
#endif

double dnorm(double x, double mu, double sigma, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x + mu + sigma;
#endif
    if(isinf(sigma)) return R_D__0;
    if(isinf(x) && mu == x) return std::numeric_limits<double>::quiet_NaN();/* x-mu is NaN */
    if (sigma <= 0) {
        if (sigma < 0){
            throw("Nan");
            return std::numeric_limits<double>::quiet_NaN();
        }
    /* sigma == 0 */
    return (x == mu) ? +INFINITY : R_D__0;
    }
    x = (x - mu) / sigma;

    if(isinf(x)) return R_D__0;

    x = fabs (x);
    if (x >= 2 * sqrt(DBL_MAX)) return R_D__0;
    if (give_log)
    return -(M_LN_SQRT_2PI + 0.5 * x * x + log(sigma));
    //  M_1_SQRT_2PI = 1 / sqrt(2 * pi)
#ifdef MATHLIB_FAST_dnorm
    // and for R <= 3.0.x and R-devel upto 2014-01-01:
    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
#else
    // more accurate, less fast :
    if (x < 5)    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;

    /* ELSE:
     * x*x  may lose upto about two digits accuracy for "large" x
     * Morten Welinder's proposal for PR#15620
     * https://bugs.r-project.org/bugzilla/show_bug.cgi?id=15620
     * -- 1 --  No hoop jumping when we underflow to zero anyway:
     *  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
     *     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
     * but "thanks" to denormalized numbers, underflow happens a bit later,
     *  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
     * for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
     * ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
     *              =IEEE=  38.58601
     * [on one x86_64 platform, effective boundary a bit lower: 38.56804]
     */
    if (x > sqrt(-2*M_LN2*(DBL_MIN_EXP + 1-DBL_MANT_DIG))) return 0.;

    /* Now, to get full accurary, split x into two parts,
     *  x = x1+x2, such that |x2| <= 2^-16.
     * Assuming that we are using IEEE doubles, that means that
     * x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).
     * If we do not have IEEE this is still an improvement over the naive formula.
     */
    double x1 = //  R_forceint(x * 65536) / 65536 =
    ldexp( R_forceint(ldexp(x, 16)), -16);
    double x2 = x - x1;
    return M_1_SQRT_2PI / sigma *
    (exp(-0.5 * x1 * x1) * exp( (-0.5 * x2 - x1) * x2 ) );
#endif
}


double attribute_hidden chebyshev_eval(double x, const double *a, const int n)
{
    double b0, b1, b2, twox;
    int i;

    if (n < 1 || n > 1000) return std::numeric_limits<double>::quiet_NaN();

    if (x < -1.1 || x > 1.1) return std::numeric_limits<double>::quiet_NaN();

    twox = x * 2;
    b2 = b1 = 0;
    b0 = 0;
    for (i = 1; i <= n; i++) {
    b2 = b1;
    b1 = b0;
    b0 = twox * b1 - b2 + a[n - i];
    }
    return (b0 - b2) * 0.5;
}


double attribute_hidden lgammacor(double x)
{
    const static double algmcs[15] = {
    +.1666389480451863247205729650822e+0,
    -.1384948176067563840732986059135e-4,
    +.9810825646924729426157171547487e-8,
    -.1809129475572494194263306266719e-10,
    +.6221098041892605227126015543416e-13,
    -.3399615005417721944303330599666e-15,
    +.2683181998482698748957538846666e-17,
    -.2868042435334643284144622399999e-19,
    +.3962837061046434803679306666666e-21,
    -.6831888753985766870111999999999e-23,
    +.1429227355942498147573333333333e-24,
    -.3547598158101070547199999999999e-26,
    +.1025680058010470912000000000000e-27,
    -.3401102254316748799999999999999e-29,
    +.1276642195630062933333333333333e-30
    };

    double tmp;

/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 *   xbig = 2 ^ 26.5
 *   xmax = DBL_MAX / 48 =  2^1020 / 3 */
#define nalgm 5
#define xbig  94906265.62425156
#define xmax  3.745194030963158e306

    if (x < 10) return std::numeric_limits<double>::quiet_NaN();
    else if (x >= xmax) {
        throw("underflow");
    /* allow to underflow below */
    }
    else if (x < xbig) {
    tmp = 10 / x;
    return chebyshev_eval(tmp * tmp * 2 - 1, algmcs, nalgm) / x;
    }
    return 1 / (x * 12);
}

double attribute_hidden stirlerr(double n)
{

#define S0 0.083333333333333333333       /* 1/12 */
#define S1 0.00277777777777777777778     /* 1/360 */
#define S2 0.00079365079365079365079365  /* 1/1260 */
#define S3 0.000595238095238095238095238 /* 1/1680 */
#define S4 0.0008417508417508417508417508/* 1/1188 */

/*
  error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
*/
    const static double sferr_halves[31] = {
    0.0, /* n=0 - wrong, place holder only */
    0.1534264097200273452913848,  /* 0.5 */
    0.0810614667953272582196702,  /* 1.0 */
    0.0548141210519176538961390,  /* 1.5 */
    0.0413406959554092940938221,  /* 2.0 */
    0.03316287351993628748511048, /* 2.5 */
    0.02767792568499833914878929, /* 3.0 */
    0.02374616365629749597132920, /* 3.5 */
    0.02079067210376509311152277, /* 4.0 */
    0.01848845053267318523077934, /* 4.5 */
    0.01664469118982119216319487, /* 5.0 */
    0.01513497322191737887351255, /* 5.5 */
    0.01387612882307074799874573, /* 6.0 */
    0.01281046524292022692424986, /* 6.5 */
    0.01189670994589177009505572, /* 7.0 */
    0.01110455975820691732662991, /* 7.5 */
    0.010411265261972096497478567, /* 8.0 */
    0.009799416126158803298389475, /* 8.5 */
    0.009255462182712732917728637, /* 9.0 */
    0.008768700134139385462952823, /* 9.5 */
    0.008330563433362871256469318, /* 10.0 */
    0.007934114564314020547248100, /* 10.5 */
    0.007573675487951840794972024, /* 11.0 */
    0.007244554301320383179543912, /* 11.5 */
    0.006942840107209529865664152, /* 12.0 */
    0.006665247032707682442354394, /* 12.5 */
    0.006408994188004207068439631, /* 13.0 */
    0.006171712263039457647532867, /* 13.5 */
    0.005951370112758847735624416, /* 14.0 */
    0.005746216513010115682023589, /* 14.5 */
    0.005554733551962801371038690  /* 15.0 */
    };
    double nn;

    if (n <= 15.0) {
    nn = n + n;
    if (nn == (int)nn) return(sferr_halves[(int)nn]);
        return(lgamma(n + 1.) - (n + 0.5)*log(n) + n - M_LN_SQRT_2PI);
    }

    nn = n*n;
    if (n>500) return((S0-S1/nn)/n);
    if (n> 80) return((S0-(S1-S2/nn)/nn)/n);
    if (n> 35) return((S0-(S1-(S2-S3/nn)/nn)/nn)/n);
    /* 15 < n <= 35 : */
    return((S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n);
}


double gammafn(double x)
{
    const static double gamcs[42] = {
    +.8571195590989331421920062399942e-2,
    +.4415381324841006757191315771652e-2,
    +.5685043681599363378632664588789e-1,
    -.4219835396418560501012500186624e-2,
    +.1326808181212460220584006796352e-2,
    -.1893024529798880432523947023886e-3,
    +.3606925327441245256578082217225e-4,
    -.6056761904460864218485548290365e-5,
    +.1055829546302283344731823509093e-5,
    -.1811967365542384048291855891166e-6,
    +.3117724964715322277790254593169e-7,
    -.5354219639019687140874081024347e-8,
    +.9193275519859588946887786825940e-9,
    -.1577941280288339761767423273953e-9,
    +.2707980622934954543266540433089e-10,
    -.4646818653825730144081661058933e-11,
    +.7973350192007419656460767175359e-12,
    -.1368078209830916025799499172309e-12,
    +.2347319486563800657233471771688e-13,
    -.4027432614949066932766570534699e-14,
    +.6910051747372100912138336975257e-15,
    -.1185584500221992907052387126192e-15,
    +.2034148542496373955201026051932e-16,
    -.3490054341717405849274012949108e-17,
    +.5987993856485305567135051066026e-18,
    -.1027378057872228074490069778431e-18,
    +.1762702816060529824942759660748e-19,
    -.3024320653735306260958772112042e-20,
    +.5188914660218397839717833550506e-21,
    -.8902770842456576692449251601066e-22,
    +.1527474068493342602274596891306e-22,
    -.2620731256187362900257328332799e-23,
    +.4496464047830538670331046570666e-24,
    -.7714712731336877911703901525333e-25,
    +.1323635453126044036486572714666e-25,
    -.2270999412942928816702313813333e-26,
    +.3896418998003991449320816639999e-27,
    -.6685198115125953327792127999999e-28,
    +.1146998663140024384347613866666e-28,
    -.1967938586345134677295103999999e-29,
    +.3376448816585338090334890666666e-30,
    -.5793070335782135784625493333333e-31
    };

    int i, n;
    double y;
    double sinpiy, value;

#ifdef NOMORE_FOR_THREADS
    static int ngam = 0;
    static double xmin = 0, xmax = 0., xsml = 0., dxrel = 0.;

    /* Initialize machine dependent constants, the first time gamma() is called.
    FIXME for threads ! */
    if (ngam == 0) {
    ngam = chebyshev_init(gamcs, 42, DBL_EPSILON/20);/*was .1*d1mach(3)*/
    gammalims(&xmin, &xmax);/*-> ./gammalims.c */
    xsml = exp(fmax2(log(DBL_MIN), -log(DBL_MAX)) + 0.01);
    /*   = exp(.01)*DBL_MIN = 2.247e-308 for IEEE */
    dxrel = sqrt(DBL_EPSILON);/*was sqrt(d1mach(4)) */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
 * (xmin, xmax) are non-trivial, see ./gammalims.c
 * xsml = exp(.01)*DBL_MIN
 * dxrel = sqrt(DBL_EPSILON) = 2 ^ -26
*/
# define ngam 22
# define xmin -170.5674972726612
# define xmax  171.61447887182298
# define xsml 2.2474362225598545e-308
# define dxrel 1.490116119384765696e-8
#endif

    if(isnan(x)) return x;

    /* If the argument is exactly zero or a negative integer
     * then return NaN. */
    if (x == 0 || (x < 0 && x == round(x))) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    y = fabs(x);

    if (y <= 10) {

    /* Compute gamma(x) for -10 <= x <= 10
     * Reduce the interval and find gamma(1 + y) for 0 <= y < 1
     * first of all. */

    n = (int) x;
    if(x < 0) --n;
    y = x - n;/* n = floor(x)  ==>    y in [ 0, 1 ) */
    --n;
    value = chebyshev_eval(y * 2 - 1, gamcs, ngam) + .9375;
    if (n == 0)
        return value;/* x = 1.dddd = 1+y */

    if (n < 0) {
        /* compute gamma(x) for -10 <= x < 1 */

        /* exact 0 or "-n" checked already above */

        /* The answer is less than half precision */
        /* because x too near a negative integer. */
        if (x < -0.5 && fabs(x - (int)(x - 0.5) / x) < dxrel) {
        throw("The answer is less than half precision because x too near a negative integer.");
        }

        /* The argument is so close to 0 that the result would overflow. */
        if (y < xsml) {
        throw ("allowed range problem, the argument is so close to 0 that the result would overflow.");
        if(x > 0) return +INFINITY;
        else return -INFINITY;
        }

        n = -n;

        for (i = 0; i < n; i++) {
        value /= (x + i);
        }
        return value;
    }
    else {
        /* gamma(x) for 2 <= x <= 10 */

        for (i = 1; i <= n; i++) {
        value *= (y + i);
        }
        return value;
    }
    }
    else {
    /* gamma(x) for     y = |x| > 10. */

    if (x > xmax) {            /* Overflow */
        throw("Overflow");
        return +INFINITY;
    }

    if (x < xmin) {            /* Underflow */
        throw("Underflow");
        return 0.;
    }

    if(y <= 50 && y == (int)y) { /* compute (n - 1)! */
        value = 1.;
        for (i = 2; i < y; i++) value *= i;
    }
    else { /* normal case */
        value = exp((y - 0.5) * log(y) - y + M_LN_SQRT_2PI +
            ((2*y == (int)2*y)? stirlerr(y) : lgammacor(y)));
    }
    if (x > 0)
        return value;

    if (fabs((x - (int)(x - 0.5))/x) < dxrel){

        /* The answer is less than half precision because */
        /* the argument is too near a negative integer. */
        throw("The answer is less than half precision because the argument is too near a negative integer.");
    }

    sinpiy = sinpi(y);
    if (sinpiy == 0) {        /* Negative integer arg - overflow */
        throw("Negative integer arg - overflow");
        return +INFINITY;
    }

    return -M_PI / (y * sinpiy * value);
    }
}



double lgammafn_sign(double x, int *sgn)
{
    double ans, y, sinpiy;

#ifdef NOMORE_FOR_THREADS
    static double xmax = 0.;
    static double dxrel = 0.;

    if (xmax == 0) {/* initialize machine dependent constants _ONCE_ */
    xmax = d1mach(2)/log(d1mach(2));/* = 2.533 e305     for IEEE double */
    dxrel = sqrt (d1mach(4));/* sqrt(Eps) ~ 1.49 e-8  for IEEE double */
    }
#else
/* For IEEE double precision DBL_EPSILON = 2^-52 = 2.220446049250313e-16 :
   xmax  = DBL_MAX / log(DBL_MAX) = 2^1024 / (1024 * log(2)) = 2^1014 / log(2)
   dxrel = sqrt(DBL_EPSILON) = 2^-26 = 5^26 * 1e-26 (is *exact* below !)
 */
#define xmax  2.5327372760800758e+305
#define dxrel 1.490116119384765625e-8
#endif

    if (sgn != NULL) *sgn = 1;

#ifdef IEEE_754
    if(ISNAN(x)) return x;
#endif

    if (sgn != NULL && x < 0 && fmod(floor(-x), 2.) == 0)
    *sgn = -1;

    if (x <= 0 && x == trunc(x)) { /* Negative integer argument */
        throw("Negative integer argument");
    return +INFINITY;/* +Inf, since lgamma(x) = log|gamma(x)| */
    }

    y = fabs(x);

    if (y < 1e-306) return -log(y); // denormalized range, R change
    if (y <= 10) return log(fabs(gammafn(x)));
    /*
      ELSE  y = |x| > 10 ---------------------- */

    if (y > xmax) {
        throw("Argument too big");
    return +INFINITY;
    }

    if (x > 0) { /* i.e. y = x > 10 */
#ifdef IEEE_754
    if(x > 1e17)
        return(x*(log(x) - 1.));
    else if(x > 4934720.)
        return(M_LN_SQRT_2PI + (x - 0.5) * log(x) - x);
    else
#endif
        return M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x);
    }
    /* else: x < -10; y = -x */
    sinpiy = fabs(sinpi(y));

    if (sinpiy == 0) { /* Negative integer argument ===
              Now UNNECESSARY: caught above */
    printf(" ** should NEVER happen! *** [lgamma.c: Neg.int, y=%g]\n",y);
        return std::numeric_limits<double>::quiet_NaN();
    }

    ans = M_LN_SQRT_PId2 + (x - 0.5) * log(y) - x - log(sinpiy) - lgammacor(y);

    if(fabs((x - trunc(x - 0.5)) * ans / x) < dxrel) {

    /* The answer is less than half precision because
     * the argument is too near a negative integer. */

    throw("The answer is less than half precision because the argument is too near a negative integer");
    }

    return ans;
}

double lgammafn(double x)
{
    return lgammafn_sign(x, NULL);
}




static double
logcf (double x, double i, double d,
       double eps /* ~ relative tolerance */)
{
    double c1 = 2 * d;
    double c2 = i + d;
    double c4 = c2 + d;
    double a1 = c2;
    double b1 = i * (c2 - i * x);
    double b2 = d * d * x;
    double a2 = c4 * c2 - b2;

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
    double c3 = c2*c2*x;
    c2 += d;
    c4 += d;
    a1 = c4 * a2 - c3 * a1;
    b1 = c4 * b2 - c3 * b1;

    c3 = c1 * c1 * x;
    c1 += d;
    c4 += d;
    a2 = c4 * a1 - c3 * a2;
    b2 = c4 * b1 - c3 * b2;

    if (fabs (b2) > scalefactor) {
        a1 /= scalefactor;
        b1 /= scalefactor;
        a2 /= scalefactor;
        b2 /= scalefactor;
    } else if (fabs (b2) < 1 / scalefactor) {
        a1 *= scalefactor;
        b1 *= scalefactor;
        a2 *= scalefactor;
        b2 *= scalefactor;
    }
    }

    return a2 / b2;
}

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
double log1pmx (double x)
{
    static const double minLog1Value = -0.79149064;

    if (x > 1 || x < minLog1Value)
    return log1p(x) - x;
    else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
        * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
        * ---------------------------------------------
        * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
       */
    double r = x / (2 + x), y = r * r;
    if (fabs(x) < 1e-2) {
        static const double two = 2;
        return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
                two / 3) * y - x);
    } else {
        static const double tol_logcf = 1e-14;
        return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
    }
    }
}


/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
double lgamma1p (double a)
{
    const double eulers_const =     0.5772156649015328606065120900824024;

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    const int N = 40;
    static const double coeffs[40] = {
    0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
    0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
    0.2058080842778454787900092413529198e-1,
    0.7385551028673985266273097291406834e-2,
    0.2890510330741523285752988298486755e-2,
    0.1192753911703260977113935692828109e-2,
    0.5096695247430424223356548135815582e-3,
    0.2231547584535793797614188036013401e-3,
    0.9945751278180853371459589003190170e-4,
    0.4492623673813314170020750240635786e-4,
    0.2050721277567069155316650397830591e-4,
    0.9439488275268395903987425104415055e-5,
    0.4374866789907487804181793223952411e-5,
    0.2039215753801366236781900709670839e-5,
    0.9551412130407419832857179772951265e-6,
    0.4492469198764566043294290331193655e-6,
    0.2120718480555466586923135901077628e-6,
    0.1004322482396809960872083050053344e-6,
    0.4769810169363980565760193417246730e-7,
    0.2271109460894316491031998116062124e-7,
    0.1083865921489695409107491757968159e-7,
    0.5183475041970046655121248647057669e-8,
    0.2483674543802478317185008663991718e-8,
    0.1192140140586091207442548202774640e-8,
    0.5731367241678862013330194857961011e-9,
    0.2759522885124233145178149692816341e-9,
    0.1330476437424448948149715720858008e-9,
    0.6422964563838100022082448087644648e-10,
    0.3104424774732227276239215783404066e-10,
    0.1502138408075414217093301048780668e-10,
    0.7275974480239079662504549924814047e-11,
    0.3527742476575915083615072228655483e-11,
    0.1711991790559617908601084114443031e-11,
    0.8315385841420284819798357793954418e-12,
    0.4042200525289440065536008957032895e-12,
    0.1966475631096616490411045679010286e-12,
    0.9573630387838555763782200936508615e-13,
    0.4664076026428374224576492565974577e-13,
    0.2273736960065972320633279596737272e-13,
    0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
    };

    const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    const double tol_logcf = 1e-14;
    double lgam;
    int i;

    if (fabs (a) >= 0.5)
    return lgammafn (a + 1);

    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
    lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);
    for (i = N - 1; i >= 0; i--)
    lgam = coeffs[i] - a * lgam;

    return (a * lgam - eulers_const) * a - log1pmx (a);
} /* lgamma1p */


double fmax2(double x, double y)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(y))
        return x + y;
#endif
    return (x < y) ? y : x;
}

/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_add (double logx, double logy)
{
    return fmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
}


/*
 * Compute the log of a difference from logs of terms, i.e.,
 *
 *     log (exp (logx) - exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
double logspace_sub (double logx, double logy)
{
    return logx + R_Log1_Exp(logy - logx);
}

/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (sum_i  exp (logx[i]) ) =
 *     log (e^M * sum_i  e^(logx[i] - M) ) =
 *     M + log( sum_i  e^(logx[i] - M)
 *
 * without causing overflows or throwing much accuracy.
 */
#ifdef HAVE_LONG_DOUBLE
# define EXP expl
# define LOG logl
#else
# define EXP exp
# define LOG log
#endif
double logspace_sum (const double* logx, int n)
{
    if(n == 0) return 0; // = log( sum(<empty>) )
    if(n == 1) return logx[0];
    if(n == 2) return logspace_add(logx[0], logx[1]);
    // else (n >= 3) :
    int i;
    // Mx := max_i log(x_i)
    double Mx = logx[0];
    for(i = 1; i < n; i++) if(Mx < logx[i]) Mx = logx[i];
    LDOUBLE s = (LDOUBLE) 0.;
    for(i = 0; i < n; i++) s += EXP(logx[i] - Mx);
    return Mx + (double) LOG(s);
}

double attribute_hidden bd0(double x, double np)
{
    double ej, s, s1, v;
    int j;

    if(isinf(x) || isinf(np) || np == 0.0) return std::numeric_limits<double>::quiet_NaN();

    if (fabs(x-np) < 0.1*(x+np)) {
    v = (x-np)/(x+np);  // might underflow to 0
    s = (x-np)*v;/* s using v -- change by MM */
    if(fabs(s) < DBL_MIN) return s;
    ej = 2*x*v;
    v = v*v;
    for (j = 1; j < 1000; j++) { /* Taylor series; 1000: no infinite loop
                    as |v| < .1,  v^2000 is "zero" */
        ej *= v;// = v^(2j+1)
        s1 = s+ej/((j<<1)+1);
        if (s1 == s) /* last term was effectively 0 */
        return s1 ;
        s = s1;
    }
    }
    /* else:  | x - np |  is not too small */
    return(x*log(x/np)+np-x);
}

double dpois_raw(double x, double lambda, int give_log)
{
    /*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
        lambda >= 0
    */
    if (lambda == 0) return( (x == 0) ? R_D__1 : R_D__0 );
    if (isinf(lambda)) return R_D__0;
    if (x < 0) return( R_D__0 );
    if (x <= lambda * DBL_MIN) return(R_D_exp(-lambda) );
    if (lambda < x * DBL_MIN) return(R_D_exp(-lambda + x*log(lambda) -lgammafn(x+1)));
    return(R_D_fexp( M_2PI*x, -stirlerr(x)-bd0(x,lambda) ));
}

/* dpois_wrap (x__1, lambda) := dpois(x__1 - 1, lambda);  where
 * dpois(k, L) := exp(-L) L^k / gamma(k+1)  {the usual Poisson probabilities}
 *
 * and  dpois*(.., give_log = TRUE) :=  log( dpois*(..) )
*/
static double
dpois_wrap (double x_plus_1, double lambda, int give_log)
{
#ifdef DEBUG_p
    REprintf (" dpois_wrap(x+1=%.14g, lambda=%.14g, log=%d)\n",
          x_plus_1, lambda, give_log);
#endif
    if (isinf(lambda))
    return 0;
    if (x_plus_1 > 1)
    return dpois_raw (x_plus_1 - 1, lambda, give_log);
    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
    return R_D_exp(-lambda - lgammafn(x_plus_1));
    else {
    double d = dpois_raw (x_plus_1, lambda, give_log);
#ifdef DEBUG_p
    REprintf ("  -> d=dpois_raw(..)=%.14g\n", d);
#endif
    return give_log
        ? d + log (x_plus_1 / lambda)
        : d * (x_plus_1 / lambda);
    }
}

/*
 * Abramowitz and Stegun 6.5.29 [right]
 */
static double
pgamma_smallx (double x, double alph, int lower_tail, int log_p)
{
    double sum = 0, c = alph, n = 0, term;

#ifdef DEBUG_p
    REprintf (" pg_smallx(x=%.12g, alph=%.12g): ", x, alph);
#endif

    /*
     * Relative to 6.5.29 all terms have been multiplied by alph
     * and the first, thus being 1, is omitted.
     */

    do {
    n++;
    c *= -x / n;
    term = c / (alph + n);
    sum += term;
    } while (fabs (term) > DBL_EPSILON * fabs (sum));

#ifdef DEBUG_p
    REprintf ("%5.0f terms --> conv.sum=%g;", n, sum);
#endif
    if (lower_tail) {
    double f1 = log_p ? log1p (sum) : 1 + sum;
    double f2;
    if (alph > 1) {
        f2 = dpois_raw (alph, x, log_p);
        f2 = log_p ? f2 + x : f2 * exp (x);
    } else if (log_p)
        f2 = alph * log (x) - lgamma1p (alph);
    else
        f2 = pow (x, alph) / exp (lgamma1p (alph));
#ifdef DEBUG_p
    REprintf (" (f1,f2)= (%g,%g)\n", f1,f2);
#endif
    return log_p ? f1 + f2 : f1 * f2;
    } else {
    double lf2 = alph * log (x) - lgamma1p (alph);
#ifdef DEBUG_p
    REprintf (" 1:%.14g  2:%.14g\n", alph * log (x), lgamma1p (alph));
    REprintf (" sum=%.14g  log(1+sum)=%.14g     lf2=%.14g\n",
          sum, log1p (sum), lf2);
#endif
    if (log_p)
        return R_Log1_Exp (log1p (sum) + lf2);
    else {
        double f1m1 = sum;
        double f2m1 = expm1 (lf2);
        return -(f1m1 + f2m1 + f1m1 * f2m1);
    }
    }
} /* pgamma_smallx() */

static double
pd_upper_series (double x, double y, int log_p)
{
    double term = x / y;
    double sum = term;

    do {
    y++;
    term *= x / y;
    sum += term;
    } while (term > sum * DBL_EPSILON);

    /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
     *       =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
     *       =  x/y * (1 + \sum_{n=1}^oo    x^n / ((y+1)*...*(y+n)))
     *       ~  x/y +  o(x/y)   {which happens when alph -> Inf}
     */
    return log_p ? log (sum) : sum;
}

/* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
static double
pd_lower_cf (double y, double d)
{
    double f= 0.0 /* -Wall */, of, f0;
    double i, c2, c3, c4,  a1, b1,  a2, b2;

#define    NEEDED_SCALE                \
      (b2 > scalefactor) {            \
        a1 /= scalefactor;            \
        b1 /= scalefactor;            \
        a2 /= scalefactor;            \
        b2 /= scalefactor;            \
    }

#define max_it 200000

#ifdef DEBUG_p
    REprintf("pd_lower_cf(y=%.14g, d=%.14g)", y, d);
#endif
    if (y == 0) return 0;

    f0 = y/d;
    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
    if(fabs(y - 1) < fabs(d) * DBL_EPSILON) { /* includes y < d = Inf */
#ifdef DEBUG_p
    REprintf(" very small 'y' -> returning (y/d)\n");
#endif
    return (f0);
    }

    if(f0 > 1.) f0 = 1.;
    c2 = y;
    c4 = d; /* original (y,d), *not* potentially scaled ones!*/

    a1 = 0; b1 = 1;
    a2 = y; b2 = d;

    while NEEDED_SCALE

    i = 0; of = -1.; /* far away */
    while (i < max_it) {

    i++;    c2--;    c3 = i * c2;    c4 += 2;
    /* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
    a1 = c4 * a2 + c3 * a1;
    b1 = c4 * b2 + c3 * b1;

    i++;    c2--;    c3 = i * c2;    c4 += 2;
    /* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
    a2 = c4 * a1 + c3 * a2;
    b2 = c4 * b1 + c3 * b2;

    if NEEDED_SCALE

    if (b2 != 0) {
        f = a2 / b2;
        /* convergence check: relative; "absolute" for very small f : */
        if (fabs (f - of) <= DBL_EPSILON * fmax2(f0, fabs(f))) {
#ifdef DEBUG_p
        REprintf(" %g iter.\n", i);
#endif
        return f;
        }
        of = f;
    }
    }

    printf(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
            f);
    return f;/* should not happen ... */
} /* pd_lower_cf() */
#undef NEEDED_SCALE


static double
pd_lower_series (double lambda, double y)
{
    double term = 1, sum = 0;

#ifdef DEBUG_p
    REprintf("pd_lower_series(lam=%.14g, y=%.14g) ...", lambda, y);
#endif
    while (y >= 1 && term > sum * DBL_EPSILON) {
    term *= y / lambda;
    sum += term;
    y--;
    }
    /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
     *       =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
     *       ~  y/lambda + o(y/lambda)
     */
#ifdef DEBUG_p
    REprintf(" done: term=%g, sum=%g, y= %g\n", term, sum, y);
#endif

    if (y != floor (y)) {
    /*
     * The series does not converge as the terms start getting
     * bigger (besides flipping sign) for y < -lambda.
     */
    double f;
#ifdef DEBUG_p
    REprintf(" y not int: add another term ");
#endif
    /* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
     *      and is unnecessary e.g. for pgamma(4e12, 121.1) */
    f = pd_lower_cf (y, lambda + 1 - y);
#ifdef DEBUG_p
    REprintf("  (= %.14g) * term = %.14g to sum %g\n", f, term * f, sum);
#endif
    sum += term * f;
    }

    return sum;
} /* pd_lower_series() */

/*
 * Compute the following ratio with higher accuracy that would be had
 * from doing it directly.
 *
 *         dnorm (x, 0, 1, FALSE)
 *       ----------------------------------
 *       pnorm (x, 0, 1, lower_tail, FALSE)
 *
 * Abramowitz & Stegun 26.2.12
 */
#define SIXTEN    16 /* Cutoff allowing exact "*" and "/" */

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)
{
/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/
    const static double a[5] = {
    2.2352520354606839287,
    161.02823106855587881,
    1067.6894854603709582,
    18154.981253343561249,
    0.065682337918207449113
    };
    const static double b[4] = {
    47.20258190468824187,
    976.09855173777669322,
    10260.932208618978205,
    45507.789335026729956
    };
    const static double c[9] = {
    0.39894151208813466764,
    8.8831497943883759412,
    93.506656132177855979,
    597.27027639480026226,
    2494.5375852903726711,
    6848.1904505362823326,
    11602.651437647350124,
    9842.7148383839780218,
    1.0765576773720192317e-8
    };
    const static double d[8] = {
    22.266688044328115691,
    235.38790178262499861,
    1519.377599407554805,
    6485.558298266760755,
    18615.571640885098091,
    34900.952721145977266,
    38912.003286093271411,
    19685.429676859990727
    };
    const static double p[6] = {
    0.21589853405795699,
    0.1274011611602473639,
    0.022235277870649807,
    0.001421619193227893466,
    2.9112874951168792e-5,
    0.02307344176494017303
    };
    const static double q[5] = {
    1.28426009614491121,
    0.468238212480865118,
    0.0659881378689285515,
    0.00378239633202758244,
    7.29751555083966205e-5
    };

    double xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
    double min = DBL_MIN;
#endif
    int i, lower, upper;

#ifdef IEEE_754
    if(ISNAN(x)) { *cum = *ccum = x; return; }
#endif

    /* Consider changing these : */
    eps = DBL_EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = fabs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
    if (y > eps) {
        xsq = x * x;
        xnum = a[4] * xsq;
        xden = xsq;
        for (i = 0; i < 3; ++i) {
        xnum = (xnum + a[i]) * xsq;
        xden = (xden + b[i]) * xsq;
        }
    } else xnum = xden = 0.0;

    temp = x * (xnum + a[3]) / (xden + b[3]);
    if(lower)  *cum = 0.5 + temp;
    if(upper) *ccum = 0.5 - temp;
    if(log_p) {
        if(lower)  *cum = log(*cum);
        if(upper) *ccum = log(*ccum);
    }
    }
    else if (y <= M_SQRT_32) {

    /* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

    xnum = c[8] * y;
    xden = y;
    for (i = 0; i < 7; ++i) {
        xnum = (xnum + c[i]) * y;
        xden = (xden + d[i]) * y;
    }
    temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)                            \
    xsq = trunc(X * SIXTEN) / SIXTEN;                \
    del = (X - xsq) * (X + xsq);                    \
    if(log_p) {                            \
        *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);    \
        if((lower && x > 0.) || (upper && x <= 0.))            \
          *ccum = log1p(-exp(-xsq * xsq * 0.5) *        \
                exp(-del * 0.5) * temp);        \
    }                                \
    else {                                \
        *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;    \
        *ccum = 1.0 - *cum;                        \
    }

#define swap_tail                        \
    if (x > 0.) {/* swap  ccum <--> cum */            \
        temp = *cum; if(lower) *cum = *ccum; *ccum = temp;    \
    }

    do_del(y);
    swap_tail;
    }

/* else      |x| > sqrt(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly     *not*    for  log_p !
 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
    else if((log_p && y < 1e170) /* avoid underflow below */
    /*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
     * Then, make use of  Abramowitz & Stegun, 26.2.13, something like
     xsq = x*x;
     if(xsq * DBL_EPSILON < 1.)
        del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
     else
        del = 0.;
     *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + log1p(-del);
     *ccum = log1p(-exp(*cum)); /.* ~ log(1) = 0 *./
      swap_tail;
     [Yes, but xsq might be infinite.]
    */
        || (lower && -37.5193 < x  &&  x < 8.2924)
        || (upper && -8.2924  < x  &&  x < 37.5193)
    ) {

    /* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
    xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
    xnum = p[5] * xsq;
    xden = xsq;
    for (i = 0; i < 4; ++i) {
        xnum = (xnum + p[i]) * xsq;
        xden = (xden + q[i]) * xsq;
    }
    temp = xsq * (xnum + p[4]) / (xden + q[4]);
    temp = (M_1_SQRT_2PI - temp) / y;

    do_del(x);
    swap_tail;
    } else { /* large x such that probs are 0 or 1 */
    if(x > 0) {    *cum = R_D__1; *ccum = R_D__0;    }
    else {            *cum = R_D__0; *ccum = R_D__1;    }
    }


#ifdef NO_DENORMS
    /* do not return "denormalized" -- we do in R */
    if(log_p) {
    if(*cum > -min)     *cum = -0.;
    if(*ccum > -min)*ccum = -0.;
    }
    else {
    if(*cum < min)     *cum = 0.;
    if(*ccum < min)    *ccum = 0.;
    }
#endif
    return;
}

static double
dpnorm (double x, int lower_tail, double lp)
{
    /*
     * So as not to repeat a pnorm call, we expect
     *
     *     lp == pnorm (x, 0, 1, lower_tail, TRUE)
     *
     * but use it only in the non-critical case where either x is small
     * or p==exp(lp) is close to 1.
     */

    if (x < 0) {
    x = -x;
    lower_tail = !lower_tail;
    }

    if (x > 10 && !lower_tail) {
    double term = 1 / x;
    double sum = term;
    double x2 = x * x;
    double i = 1;

    do {
        term *= -i / x2;
        sum += term;
        i += 2;
    } while (fabs (term) > DBL_EPSILON * sum);

    return 1 / sum;
    } else {
    double d = dnorm (x, 0., 1., false);
    return d / exp (lp);
    }
}

/*
 * Asymptotic expansion to calculate the probability that Poisson variate
 * has value <= x.
 * Various assertions about this are made (without proof) at
 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */

double pnorm(double x, double mu, double sigma, int lower_tail, int log_p)
{
    double p, cp;

    /* Note: The structure of these checks has been carefully thought through.
     * For example, if x == mu and sigma == 0, we get the correct answer 1.
     */
#ifdef IEEE_754
    if(ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
    return x + mu + sigma;
#endif
    if(isinf(x) && mu == x) return std::numeric_limits<double>::quiet_NaN();/* x-mu is NaN */
    if (sigma <= 0) {
        if(sigma < 0){
            throw("error sigma < 0");
            return std::numeric_limits<double>::quiet_NaN();
        };
    /* sigma = 0 : */
    return (x < mu) ? R_DT_0 : R_DT_1;
    }
    p = (x - mu) / sigma;
    if(isinf(p))
    return (x < mu) ? R_DT_0 : R_DT_1;
    x = p;

    pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

    return(lower_tail ? p : cp);
}


static double
ppois_asymp (double x, double lambda, int lower_tail, int log_p)
{
    static const double coefs_a[8] = {
    -1e99, /* placeholder used for 1-indexing */
    2/3.,
    -4/135.,
    8/2835.,
    16/8505.,
    -8992/12629925.,
    -334144/492567075.,
    698752/1477701225.
    };

    static const double coefs_b[8] = {
    -1e99, /* placeholder */
    1/12.,
    1/288.,
    -139/51840.,
    -571/2488320.,
    163879/209018880.,
    5246819/75246796800.,
    -534703531/902961561600.
    };

    double elfb, elfb_term;
    double res12, res1_term, res1_ig, res2_term, res2_ig;
    double dfm, pt_, s2pt, f, np;
    int i;

    dfm = lambda - x;
    /* If lambda is large, the distribution is highly concentrated
       about lambda.  So representation error in x or lambda can lead
       to arbitrarily large values of pt_ and hence divergence of the
       coefficients of this approximation.
    */
    pt_ = - log1pmx (dfm / x);
    s2pt = sqrt (2 * x * pt_);
    if (dfm < 0) s2pt = -s2pt;

    res12 = 0;
    res1_ig = res1_term = sqrt (x);
    res2_ig = res2_term = s2pt;
    for (i = 1; i < 8; i++) {
    res12 += res1_ig * coefs_a[i];
    res12 += res2_ig * coefs_b[i];
    res1_term *= pt_ / i ;
    res2_term *= 2 * pt_ / (2 * i + 1);
    res1_ig = res1_ig / x + res1_term;
    res2_ig = res2_ig / x + res2_term;
    }

    elfb = x;
    elfb_term = 1;
    for (i = 1; i < 8; i++) {
    elfb += elfb_term * coefs_b[i];
    elfb_term /= x;
    }
    if (!lower_tail) elfb = -elfb;
#ifdef DEBUG_p
    REprintf ("res12 = %.14g   elfb=%.14g\n", elfb, res12);
#endif

    f = res12 / elfb;

    np = pnorm (s2pt, 0.0, 1.0, !lower_tail, log_p);

    if (log_p) {
    double n_d_over_p = dpnorm (s2pt, !lower_tail, np);
#ifdef DEBUG_p
    REprintf ("pp*_asymp(): f=%.14g     np=e^%.14g  nd/np=%.14g  f*nd/np=%.14g\n",
          f, np, n_d_over_p, f * n_d_over_p);
#endif
    return np + log1p (f * n_d_over_p);
    } else {
    double nd = dnorm (s2pt, 0., 1., log_p);

#ifdef DEBUG_p
    REprintf ("pp*_asymp(): f=%.14g     np=%.14g  nd=%.14g  f*nd=%.14g\n",
          f, np, nd, f * nd);
#endif
    return np + f * nd;
    }
} /* ppois_asymp() */


double pgamma_raw (double x, double alph, int lower_tail, int log_p)
{
/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

    double res;

#ifdef DEBUG_p
    REprintf("pgamma_raw(x=%.14g, alph=%.14g, low=%d, log=%d)\n",
         x, alph, lower_tail, log_p);
#endif
    R_P_bounds_01(x, 0., ML_POSINF);

    if (x < 1) {
    res = pgamma_smallx (x, alph, lower_tail, log_p);
    } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
    /* incl. large alph compared to x */
    double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
    double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
    REprintf(" alph 'large': sum=pd_upper*()= %.12g, d=dpois_w(*)= %.12g\n",
         sum, d);
#endif
    if (!lower_tail)
        res = log_p
        ? R_Log1_Exp (d + sum)
        : 1 - d * sum;
    else
        res = log_p ? sum + d : sum * d;
    } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
    /* incl. large x compared to alph */
    double sum;
    double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
    REprintf(" x 'large': d=dpois_w(*)= %.14g ", d);
#endif
    if (alph < 1) {
        if (x * DBL_EPSILON > 1 - alph)
        sum = R_D__1;
        else {
        double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
        /* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
        sum = log_p ? log (f) : f;
        }
    } else {
        sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
        sum = log_p ? log1p (sum) : 1 + sum;
    }
#ifdef DEBUG_p
    REprintf(", sum= %.14g\n", sum);
#endif
    if (!lower_tail)
        res = log_p ? sum + d : sum * d;
    else
        res = log_p
        ? R_Log1_Exp (d + sum)
        : 1 - d * sum;
    } else { /* x >= 1 and x fairly near alph. */
#ifdef DEBUG_p
    REprintf(" using ppois_asymp()\n");
#endif
    res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
    }

    /*
     * We lose a fair amount of accuracy to underflow in the cases
     * where the final result is very close to DBL_MIN.     In those
     * cases, simply redo via log space.
     */
    if (!log_p && res < DBL_MIN / DBL_EPSILON) {
    /* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
#ifdef DEBUG_p
    REprintf(" very small res=%.14g; -> recompute via log\n", res);
#endif
    return exp (pgamma_raw (x, alph, lower_tail, 1));
    } else
    return res;
}


double pgamma(double x, double alph, double scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(alph) || ISNAN(scale))
    return x + alph + scale;
#endif
    if(alph < 0. || scale <= 0.){
        throw("Error bonk!");
        return std::numeric_limits<double>::quiet_NaN();
    }
    x /= scale;
#ifdef IEEE_754
    if (ISNAN(x)) /* eg. original x = scale = +Inf */
    return x;
#endif
    if(alph == 0.) /* limit case; useful e.g. in pnchisq() */
    return (x <= 0) ? R_DT_0: R_DT_1; /* <= assert  pgamma(0,0) ==> 0 */
    return pgamma_raw (x, alph, lower_tail, log_p);
}
/* From: terra@gnome.org (Morten Welinder)
 * To: R-bugs@biostat.ku.dk
 * Cc: maechler@stat.math.ethz.ch
 * Subject: Re: [Rd] pgamma discontinuity (PR#7307)
 * Date: Tue, 11 Jan 2005 13:57:26 -0500 (EST)

 * this version of pgamma appears to be quite good and certainly a vast
 * improvement over current R code.  (I last looked at 2.0.1)  Apart from
 * type naming, this is what I have been using for Gnumeric 1.4.1.

 * This could be included into R as-is, but you might want to benefit from
 * making logcf, log1pmx, lgamma1p, and possibly logspace_add/logspace_sub
 * available to other parts of R.

 * MM: I've not (yet?) taken  logcf(), but the other four
 */

/*
int main(int argc, const char * argv[]) {
    double x=10;
    double df=135;
    int lower_tail=0;
    int log_p=1;
    std::cout << pgamma(x, df/2., 2., lower_tail, log_p) << "\n";
}
*/
