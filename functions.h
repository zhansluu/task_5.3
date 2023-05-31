#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED
#include <cmath>

double f(double x)
{
    return sin(x);
}

double rho(double x)
{
    return pow(x, -0.5);
}

double phi(double x)
{
    return sin(x)*pow(x, -0.5);
}
// находит коэффициенты многочленов Лежандра
void fill_legendre(int degree, double** polynom_coeff)
{
    polynom_coeff[0][0] = 1;
    polynom_coeff[1][1] = 1;
    for (int i = 2; i <= degree; i++) {
        for (int j = 1; j <= degree; j++) {
            polynom_coeff[i][j] = (2 * i - 1) * polynom_coeff[i - 1][j - 1] / i - (i - 1) * polynom_coeff[i - 2][j] / i;
        }
        polynom_coeff[i][0] = polynom_coeff[i - 2][0] * (1 - i) / i;
    }
}

double legendre(double x, size_t N)
{
    double P;
    if (N == 0) P = 1;
    else if (N == 1) P = x;
    else if (N > 1)
        P = ((2*N-1)*legendre(x, N-1)*x)/N - ((N-1)*legendre(x, N-2))/N;
    return P;
}
// вычисляет многочлен лежандра степени degree в точке x
/*//double legendre(double x, int degree, double** polynom_coeff) {
    if (degree == 0) {
        return 1;
    }
    else if (degree == 1) {
        return x;
    }
    else
    {
        double value = 0;
        for (int i = 0; i <= degree; i++)
        {
            value += polynom_coeff[degree][i] * pow(x, i);
        }
        return value;
    }
}*/

/*double legendre(double x, size_t N)
{
    double P;
    if (N == 0) P = 1;
    else if (N == 1) P = x;
    else if (N > 1)
        P = ((2*N-1)*legendre(x, N-1)*x)/N - ((N-1)*legendre(x, N-2))/N;
    return P;
}*/

/*Метод бисекции (половинного деления)

    Делим отрезок пополам, пока его длина не будет меньше 2*epsilon.
    Как только мы достигли этого условия, берем X, который находится на середине получившегося отрезка.*/
    double bisection_method (double a, double b, double epsilon, size_t N)
    {
        double c, X/*, delta*/;
        size_t counter = 0;

        while ((b - a) > 2*epsilon)
        {
            c = (a+b)/2;
            if (legendre(a, N)*legendre(c, N) <= 0)
                b = c;
            else
                a = c;
            counter++;
        }
        X = (a + b)/2;
        //delta = (b - a)/2;
        return X;
        /*std::cout << "x = " << std::setprecision(15) << X << std::endl;
        std::cout << "Длина последнего отрезка: " << delta << std::endl;
        std::cout << "Абсолютная величина невязки для приближенного решения: " << fabs(legendre(X) - 0) << std::endl;
        std::cout << "Количество шагов для достижения точности эпсилон: " << counter << std::endl << std::endl;*/
    }

/*Метод бисекции (половинного деления)

    Делим отрезок пополам, пока его длина не будет меньше 2*epsilon.
    Как только мы достигли этого условия, берем X, который находится на середине получившегося отрезка.*/
 /*   double bisection_method (double a, double b, double epsilon, size_t N, double** polynom_coeff)
    {
        double c, X;
        size_t counter = 0;

        while ((b - a) > 2*epsilon)
        {
            c = (a+b)/2;
            if (legendre(a, N, polynom_coeff)*legendre(c, N, polynom_coeff) <= 0)
                b = c;
            else
                a = c;
            counter++;
        }
        X = (a + b)/2;
        return X;

    }
*/
#endif // FUNCTIONS_H_INCLUDED
