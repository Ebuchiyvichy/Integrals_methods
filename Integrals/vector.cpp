#include "MatrixClass.h"

void print(std::vector<double> x)
{
	for (int i = 0; i != x.size(); i++)
		std::cout << x[i] << '\t';
	std::cout << std::endl;
}
void print_in_file(std::vector<double> x, std::ofstream fout)
{
	for (int i = 0; i != x.size(); i++)
		fout << x[i] << '\t';
	fout << std::endl;
}

std::vector<double>	cpy_vector(std::vector<double> tmp, std::vector<double> x)
{

	for (int i = 0; i < x.size(); i++)
		tmp[i] = x[i];
	return (tmp);
}

// ���������������� �������� ��� �������
std::vector<double> operator * (double a, std::vector<double> b)
{
	std::vector<double> c(b);
	for (int i = 0; i != c.size(); i++)
		c[i] = a * b[i];
	return c;
}
std::vector<double> operator + (std::vector<double> a, std::vector<double> b)
{
	std::vector<double> c(b);
	for (int i = 0; i != c.size(); i++)
		c[i] = a[i] + b[i];
	return c;
}
std::vector<double> operator - (std::vector<double> a, std::vector<double> b)
{
	std::vector<double> c(b);
	for (int i = 0; i != c.size(); i++)
		c[i] = a[i] - b[i];
	return c;
}
std::vector<double> operator / (std::vector<double> a, double b)
{
	std::vector<double> c(a);
	for (int i = 0; i != c.size(); i++)
		c[i] = a[i] / b;
	return c;
}

// ����� ��������
double norm(std::vector<double> x)
{
	double sum = 0;
	for (int i = 0; i != x.size(); i++)
		sum += (x[i] * x[i]);
	return sqrt(sum);
}
double	norm(std::vector<double> a, std::vector<double> b)
{
	double	max;
	max = fabs(a[0] - b[0]);
	for (int i = 1; i < a.size(); i++)
		if (fabs(a[i] - b[i]) > max)
			max = fabs(a[i] - b[i]);
	return max;
}