#include <iostream>
#include "functions.h"
#include <iomanip>

using namespace std;

int main()
{
    setlocale(0, "Russian");
    cout << "������� 5.3. ������������ ���������� ��������� ��� ������ ��������� �� ������" << endl << endl;

    double a, b, a1 = -1, b1 = 1;
    size_t N, m;
    cout << "������� ������ ������ ��������������: � = ";
    cin >> a;
    cout << "������� ������� ������ ��������������: b = ";
    cin >> b;
    cout << "������� ���������� �����: N = ";
    cin >> N;
    while (N <= 0)
    {
        cout << "������������ �������� N! ������� ����� ��������:";
        cin >> N;
    }
    cout << "������� ����� ����������� ������� [a,b]: m = ";
    cin >> m;
    while (m <= 0)
    {
        cout << "������������ �������� m! ������� ����� ��������:";
        cin >> m;
    }
    double h, x1[10000], x2[10000], y1[10000], y2[10000], x[100], coeff[100];
        size_t counter = 0;
        h = (b1-a1)/10000;
        x1[counter] = a1;
        x2[counter] = x1[counter] + h;
        while (x2[counter] <= b1)
        {
            y1[counter] = legendre(x1[counter], N);
            y2[counter] = legendre(x2[counter], N);
            if (y1[counter]*y2[counter] <= 0)
            {
                counter++;
                x1[counter] = x2[counter-1];
                x2[counter] = x2[counter-1] + h;
            }
            else
            {
                x1[counter] = x2[counter];
                x2[counter] = x2[counter] + h;
            }
        }
        cout << "���� �������� �� ������:" << endl;
        //cout << "����� �������� �������� �����: " << counter << endl << endl;
        size_t counter_2 = 0, i = 1;
        do
        {
            x[i] = bisection_method(x1[counter_2], x2[counter_2], 0.000000000001, N);
            cout << "x[" << i << "] = " << x[i] << "; ";
            counter_2 ++;
            i++;
        } while (counter_2 != counter);

        //������������ �� ������
        for(size_t k = 1; k <= N; k++)
        {
            coeff[k] = (2*(1-x[k]*x[k]))/(N*N*legendre(x[k],N-1)*legendre(x[k],N-1));
            if (coeff[k] <= 0) { cout << "����������� �����������!"; return 0;}
        }

        cout << endl;
        double sum = 0;
        h = (b-a)/m;
        for (size_t j = 0; j < m; j++)
        {
            for (size_t k = 1; k <= N; k++)
            {
                sum += coeff[k] * phi(h * x[k] / 2 + a + j * h + h / 2);
            }
        }

        sum *= (h / 2);

        cout << "���������� �������� ���������: " << setprecision(12) << sum;

    return 0;
}
