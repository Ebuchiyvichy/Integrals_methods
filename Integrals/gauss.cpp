
#include "MatrixClass.h"
void delim(Matrix  &A, int k, std::vector<double>& x)
{
	int i = k;
	double temp = A.value[i][k];
	for (int j = k; j != A.size; j++)
	{
		A.value[i][j] = A.value[i][j] / temp;
	}
	x[i] = x[i] / temp;
}

void vych(Matrix &A, int k, std::vector<double>& x)
{
	for (int i = (k + 1); i != A.size; i++)
	{
		//columns
		double temp = A.value[i][k];
		for (int j = 0; j != A.size; j++)
			A.value[i][j] = A.value[i][j] - temp * A.value[k][j];
		x[i] -= temp * x[k];
	}
}

void GaussRight(Matrix &A, std::vector<double>& x)
{
	// for k-string our matrix
	for (int k = 0; k != A.size; k++)
	{
		//search max in column
		double max = A.value[k][k];
		int maxstring = k;
		for (int i = k; i != A.size; i++)
		{
			if (A.value[i][k] > max)
			{
				max = A.value[i][k];
				maxstring = i;
			}
		}
		//change rvalue
		double temp = x[k];
		x[k] = x[maxstring];
		x[maxstring] = temp;
		// change strings
		std::swap(A.value[k], A.value[maxstring]);
		delim(A, k, x);
		vych(A, k, x);
	}
	//std::cout << "Matrix A after straight run of Gauss:" << std::endl;
	//A.print();
}

void GaussLeft(const Matrix &A, std::vector <double>& x)
{
	for (int i = A.size - 1; i >= 0; i--)
	{
		for (int j = i + 1; j < A.size; j++)
			x[i] = x[i] - A.value[i][j] * x[j];
		x[i] = x[i] / A.value[i][i];
	}
	// std::cout << "Vector in Gauss x:" << std::endl;
	// for (int i = 0; i != A.size; i++)
	//   std::cout << std::setw(8) << x[i] << std::endl;
	// std::cout << std::endl;
}

double	norm(std::vector <double> x, int DIM)
{
	double	norm = 0;

	for (int i = 0; i < DIM; i++)
		norm += x[i] * x[i];
	return (sqrt(norm));
}

void hod(const unsigned int DIM, std::vector<double>& b, Matrix & a) //�������� ��� ������
{
	for (int i = (DIM - 1); i >= 0; i--)
	{
		for (int j = i + 1; j < DIM; j++)
		{
			b[i] = b[i] - a.value[i][j] * b[j]; //��������� �� ������� ������ ����� ���� ��������� ������� ���� �� ������ ����� �������� �� �������� ���������
		}
		b[i] = b[i] / a.value[i][i]; //������� �������� ������� ������ ����� �� �������� ������� ���������
	}
}

bool findx(const unsigned int DIM, std::vector<double>& b, Matrix& a) //������ ��� ������ � ��������� �� �������������
{
	double eps;     //��������������� ���������� ��� ��������� ��������� (���������� max)
	double* vspom;  //��������������� ������ ��� ������������ ����� �������
	double vspom2;   //��������������� ���������� ��� ������������ ��������� �������
	bool flag = true;  //���������� ����������(�������� ������������� �������)
	double EPS = 1e-20;

	for (int k = 0, h; k < (DIM); k++)
	{
		eps = fabs(a.value[k][k]);  //������ ������������� �������� � ���������� ���������
									//������ ������ �� ������ ������������� ��������(���� ������� � ����� �� ���� �� ����������)
		h = k;   //������ ������ ������ (�� �������, ��������� ����)

		if (fabs(a.value[k][k]) < EPS) { flag = true; } //���� ������ ������� �������, �� ������������ ����
		else { flag = false; } //���� �� �������, �� �������������� ����

		for (int i = (k + 1); i < DIM; i++)
		{
			if (fabs(a.value[i][k]) > eps && fabs(a.value[i][k]) > EPS) //���� ������� ������ ����������� max � ������� � �� ����
			{
				eps = a.value[i][k];  h = i; flag = false;
			}; //��������� ������, �������� max � �������������� ����
		}

		if (!flag) //���� ������� �� �������
		{
			std::swap(a.value[h], a.value[k]);  //����� ��������
			std::swap(b[h], b[k]);//����� ���������� �������
		}
		else { std::cout << "������� ���������\n";   return flag; }

		for (int i = k + 1; i < DIM; i++)
		{
			b[i] = b[i] - b[k] / a.value[k][k] * a.value[i][k]; //��������� �� ��������� �������
			for (int j = (DIM - 1); j >= 0; j--)
			{
				a.value[i][j] = a.value[i][j] - a.value[k][j] / a.value[k][k] * a.value[i][k]; //��������� ���� ��������� ��� ������������
			}
		}

	}
	return flag;
}