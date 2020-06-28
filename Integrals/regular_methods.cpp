#include "MatrixClass.h"

//метод квадратур с трапециями
double	qudrature_method(double N, double a, double b)
{
	double	h = (b - a) / (N - 1);
	Matrix	A(N);
	std::vector<double>	r(N);
	std::ofstream fout;
	double	error = 0;

	fout.open("quadrature.txt");

	for (int i = 0; i < N; i++)
	{
		A.value[i][0] = -h / 2 * K(a + i * h, a);
		for (int k = 1; k < N - 1; k++)
			A.value[i][k] = -h * K(a + i * h, a + k * h);
		A.value[i][N - 1] = -h / 2 * K(a + i * h, b);
		A.value[i][i] += 1;
	}

	for (int i = 0; i < N; i++)
		r[i] = f(a + h * i);

	findx(N, r, A);
	hod(N, r, A);
	for (int i = 0; i < N; i++)
		fout << a + i * h << '\t' << r[i] << '\n';

	for (int i = 0; i < N; i++)	//ошибка для первого теста
		if (fabs(1.0 - r[i]) > error)
			error = fabs(1.0 - r[i]);

	fout.close();
	return error;
}
//метод простой итерации с трапециями
void	simple_iteration_method(double N, double a, double b)
{
	double				h = (b - a) / (N - 1);
	std::vector<double>	u0(N);
	std::vector<double>	uk(N);
	std::vector<double>	func(N);
	std::vector<double>	Kore(N);
	std::ofstream		fout;
	double				error = 0;


	int					iter = 0;

	fout.open("simple_iteration.txt");

	for (int i = 0; i < N; i++)
		func[i] = f(a + i * h);
	for (int i = 0; i < N; i++)
		uk[i] = 0;

	do {
		iter++;
		u0 = uk;
		for (int i = 0; i < N; i++) {
			Kore[i] = 0;
			for (int k = 0; k < N - 1; k++)
				Kore[i] += h / 2 * (K(a + i * h, a + (k) * h) * u0[k] + K(a + i * h, a + (k+1) * h) * u0[k+1]);
		}
		uk = func + Kore;
	} while (norm(u0, uk) > EPS);

	std::cout << iter << " iterations to " << EPS << " precision\n";

	for (int i = 0; i < N; i++)
		fout << i * h << '\t' << uk[i] << '\n';
	fout.close();
	for (int i = 0; i < N; i++)	//ошибка для первого теста
		if (fabs(1.0 - uk[i]) > error)
			error = fabs(1.0 - uk[i]);
	std::cout << "Error is " << error << '\n';
}
//вычисление интеграла для вырожденного ядра
double	integral(int N, double a, double b, int i, int j)
{
	double I;
	double	h = (b - a) / (N - 1);

	I = 0;
	for (int k = 0; k < N - 1; k++)
		I += h / 2 * (psi(a + h * (k), i) * phi(a + h * (k), j) + psi(a + h * (k+1), i) * phi(a + (k+1) * h, j));

	return I;
}
//метод при вырожденном ядре
void	degenerate_core(double N, double a, double b, int m)
{
	double	h = (b - a) / (N - 1);
	std::vector<double>	u(N);
	std::vector<double>	C(m);
	Matrix				A(m);
	std::ofstream		fout;
	double				error = 0;

	fout.open("degenerate_core.txt");

	
	//ищем betta
	for (int i = 0; i < m; i++) {
		C[i] = 0;
		for (int k = 0; k < N - 1; k++)
			C[i] += h / 2 * (psi(a + h * (k), i) * f(a + h * (k)) + psi(a + h * (k+1), i) * f(a + h * (k+1)));
	}
	//составляем матрицу
	for (int i = 0; i < m; i++) {
		for (int k = 0; k < m; k++)
			A.value[i][k] = -integral(N, a, b, i, k);
		A.value[i][i] += 1;
	}
	findx(m, C, A);
	hod(m, C, A);

	for (int k = 0; k < N; k++) {
		u[k] = f(a + h * k);
		for (int i = 0; i < m; i++)
			u[k] += C[i] * phi(a + k * h, i);
	}
	for (int i = 0; i < N; i++)
		fout << i * h << '\t' << u[i] << '\n';
	fout.close();

	for (int i = 0; i < N; i++)	//ошибка для первого теста
		if (fabs(1.0 - u[i]) > error)
			error = fabs(1.0 - u[i]);
	std::cout << "Error is " << error << '\n';
}

void	singular_core(double N)
{
	double				l = 2 * PI / N;
	std::vector<double>	ki(N);
	std::vector<double>	ci(N);
	std::vector<double>	kj(N);
	std::vector<double>	cj(N);
	std::vector<double>	b(N + 1);
	Matrix				A(N + 1);
	std::ofstream		fout;
	double				error = 0;

	fout.open("singular.txt");

	for (int i = 0; i < N; i++)
	{
		ki[i] = cos(l * (i + 0.5));	kj[i] = sin(l * (i + 0.5));
		ci[i] = cos(l * i);			cj[i] = sin(l * i);
	}
	for (int i = 0; i < N; i++)
		b[i] = f_singular((i + 0.5) * l);

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
			A.value[i][j] = (kj[i] * (ki[i] - ci[j]) - ki[i] * (kj[i] - cj[j])) / N / (pow(ki[i] - ci[j], 2) + pow(kj[i] - cj[j], 2));
		A.value[i][N] = 1.;	A.value[N][i] = 1.;
	}
	A.value[N][N] = 0.;	b[N] = 0.;

	findx(N + 1, b, A);
	hod(N + 1, b, A);

	for (int i = 0; i < N; i++)
		fout << ci[i] << '\t' << cj[i] << '\t' << b[i] << '\n';
	fout.close();

	for (int i = 0; i < N - 1; i++)
		/*if (fabs(b[i] - 2*sin(i*l)) > error)//variant2
			error = fabs(b[i] - 2*sin(i*l));*/
		if (fabs(b[i] - 2 * sin(3*i*l)) > error)//variant6
			error = fabs(b[i] - 2 * sin(3*i*l));
	std::cout << error << '\n';

	std::cout << "R is " << b[N] << '\n';
}