/* Вычисление корней системы нелинейных урнавений методом Ньютона */

// Для функций sin и cos
#include <math.h>
// Для printf
#include <stdio.h>
// malloc
#include <stdlib.h>

// Размер системы - сколько x и соответственно сколько уравнений
#define NSIZE 100

// Начальное приближение x(0)
#define XSTART 0.87
// Забегая вперед это начальное приближение дает ответ х=0
// Если задать начальное приближение больше чем 2Пи то ответ будет 2Пи и тд(тк ф-я периодическая)

// Уровень малости невязки (мин норма вектора y)
#define EPS 1e-9

// В общем виде систему представим как
// func(x, i) = 0, где x - вектор x от 1 до NSIZE, i - номер уравнения в системе(от 0)
double	func(double *x, int i)
{
	return (cos(x[i]) - 1);
}

// Производная по x(j), где i - номер уравнения, j - номер переменной(от 0)
double	dfunc(double *x, int i, int j)
{
	if (i == j)
		return (-sin(x[j]));
	else
		return (0);
}

// Системная функция очитик памяти матрицы
void	matrix_free(double **matrix)
{
	int		i;

	i = -1;
	while (++i < NSIZE)
		if (matrix[i])
			free(matrix[i]);
	free(matrix);
}

// Создание и инициализаия вектора
double	*vec_init(double x0)
{
	double	*res;
	int		i;

	i = -1;
	res = (double*)malloc(sizeof(double) * NSIZE);
	if (!res)
		return (NULL);
	while (++i < NSIZE)
		res[i] = x0;
	return (res);
}

// Норма вектора
double	vec_norm(double *y)
{
	int		i;
	double	sum;
	
	sum = 0;
	i = -1;
	while (++i < NSIZE)
		sum += y[i] * y[i];
	return (sqrt(sum));
}

// Создание матрицы
double	**matrix_create(void)
{
	int		i;
	double	**matrix;

	i = -1;
	matrix = (double**)malloc(sizeof(double*) * NSIZE);
	if (!matrix)
		return (NULL);
	while (++i < NSIZE)
	{
		matrix[i] = (double*)malloc(sizeof(double) * NSIZE);
		if (!matrix[i])
			return (NULL);
	}
	return (matrix);
}

// Обнуление матрицы
void	null_matrix(double **matrix)
{
	int		i;
	int		j;

	i = -1;
	while (++i < NSIZE)
	{
		j = -1;
		while (++j < NSIZE)
			matrix[i][j] = 0;
	}
}

// Нахождение обратной диагональной матрицы (на месте)
void	obr_diag_matrix(double **matrix)
{
	int		i;
	int		j;
	int		flag;

	i = -1;
	flag = 0; // Для контроля диагональности
	while (++i < NSIZE)
	{
		j = -1;
		while (++j < NSIZE)
			if (i == j)
				matrix[i][j] = 1 / matrix[i][j];
			else if (matrix[i][j]) // Проверяем точно ли матрица диагональнаяш
				flag = 1;
	}
	if (flag)
	{
		printf("Ошибка матрицы\n");
		null_matrix(matrix);
	}
}

// Умножение матрицы на вектор (результат в вектор res)
void	matrix_on_vec(double **matrix, double *vec, double *res)
{
	int		i;
	int		j;

	i = -1;
	while (++i < NSIZE)
	{
		j = -1;
		while (++j < NSIZE)
			res[i] += matrix[i][j] * vec[j];
	}
}

// Отрицательный вектор(на месте)
void	neg_vec(double *vec)
{
	int		i;

	i = -1;
	while (++i < NSIZE)
		vec[i] = -vec[i];
}

// Приращение вектора х (на месте)
void	add_vec(double *x, double *delta)
{
	int		i;

	i = -1;
	while (++i < NSIZE)
		x[i] += delta[i];
}

// Вычиление вектора у по известным х, результат во второй  аргумент res
void	vec_func(double *x, double *res)
{
	int		i;

	i = -1;
	while (++i < NSIZE)
		res[i] = func(x, i);
}

// Вычисление матрицы Якобиан, результат во второй аргумент matrix
void	jakobian(double *x, double **matrix)
{
	int		i;
	int		j;

	i = -1;
	while (++i < NSIZE)
	{
		j = -1;
		while (++j < NSIZE)
			matrix[i][j] = dfunc(x, i, j);
	}
}

// Печать матрицы (длдя отладки)
void	print_matrix(double **matrix)
{
	int		i;
	int		j;

	i = -1;
	while (++i < NSIZE)
	{
		j = -1;
		while (++j < NSIZE)
			printf("%f ", matrix[i][j]);
		printf("\n");
	}
}

// Печать вектора (для отладки)
void	print_vec(double *vec)
{
	int		i;

	i = -1;
	while (++i < NSIZE)
		printf("%f ", vec[i]);
	printf("\n");
}

int	main(void)
{
	double	*x; //  Вектор х
	double	*y; // Вектор у, должен получиться нулевым
	double	*delta; // Вектор приращений х
	double	**jmatrix; // Матрица Якоби (до и после обращения)
	int		i; // Используем итератор для контроля кол-ва итераций
	double	y_norm;

	// Создаем и инициализируем вектор х начальным приближением
	x = vec_init(XSTART);
	// Создаем искомый вектор приращений х
	delta = vec_init(0);
	// Создаем рабочую матрицу (чанала в ней будет матрица Якоби, потом она же обращенная)
	jmatrix = matrix_create();
	// Создаем рабочий вектор у, значения в нем будут меняться на каждой итерации
	y = vec_init(0);
	if (!x || !delta || !jmatrix || !y)
		return (1);
	i = -1;
	// Вычисляем значение вектора у для начальных приближений
	vec_func(x, y);
	y_norm = vec_norm(y);
	while (y_norm >= EPS)
	{
		//printf("Iteration: %d\n", i); // Отладка
		//print_vec(y);// Отладка
		// Вычисляем якобиан
		null_matrix(jmatrix);
		jakobian(x, jmatrix);
		//print_matrix(jmatrix); // Отладка 
		// Обратный якобиан
		obr_diag_matrix(jmatrix);
		//print_matrix(jmatrix); // Отладка 
		// Отрицание вектора у
		neg_vec(y);
		// Получаем вектор приращений х путем умножения отрицательного вектора у на обратный якобиан
		matrix_on_vec(jmatrix, y, delta);
		// Вычисляем новый х
		add_vec(x, delta);
		//print_vec(x); // Отладка
		// Вычисляем новый у
		vec_func(x, y);
		//print_vec(y);// Отладка
		y_norm = vec_norm(y);
		i++;
	}
	printf("Итераций: %d\n", i);
	printf("Норма вектора у: %.10f\n", y_norm);
	i = -1;
	while (++i <  NSIZE)
		printf("x[%d] = %f\n", i, x[i]);
	free(x);
	free(delta);
	free(y);
	matrix_free(jmatrix);
	return (0);
}
