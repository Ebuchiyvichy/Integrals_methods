#include "MatrixClass.h"

double	EPS = 1.e-6;
double	PI = 3.1415926535;


double	K(double x, double s)
{
//	return 1. / 2 * (1. - x * cos(x * s));//tests
//	return (1. - x * cos(x * s));//variant2
	return x * (cos(x*s) - exp(-x * s));//variant6
}
double	f(double x)
{
//	return 1. / 2 * (1 + sin(x));//test1
	return x * x + sqrt(x);//test2
//	return sin(x);//variant2test1
//	return 2 - exp(-x) - sin(x);//variant6test1
}
double	phi(double x, int i)
{
	//test1
	if (i == 0)
		return 1. / 2 * (1 - x);
	else if (i == 1)
		return 1. / 4 * pow(x, 3);
	else if (i == 2)
		return -1. / 45 * pow(x, 5);
	else if (i == 3)
		return 1. / 1440 * pow(x, 7);
}
double	psi(double x, int i)
{
	//test1
	return pow(x, i * 2);
}
double	f_singular(double phi)
{
//	return cos(phi);//variant2
	return cos(3 * phi);//variant6
}

int	main()
{
	//пределы интегрирования
//	double	a = 0.;
	double	a = 0.1;
	double	b = 1.;
	int		N = 100;
	std::ofstream	fout;

	qudrature_method(100, a, b);
	std::cout << "quadrature method is over\n";

//	simple_iteration_method(100, a, b);
//	std::cout << "simple iteration method is over\n";

//	degenerate_core(100, a, b, 4);
//	std::cout << "method for degenerate core is over\n";

//	singular_core(10);
//	std::cout << "method with singular core is over\n";

//	fout.open("error.txt");
//	for (int i = 0; i < 5; i++) {
//		fout << qudrature_method(N * pow(2, i), a, b) << '\n';
//		std::cout << i << "is over\n";
//	}

//	fout.close();
	system("pause");

	return 0;
}