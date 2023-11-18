#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

using namespace std;

const double g = 9.81;



/// @brief Физические свойства трубы
struct simple_pipe_t
{
	/// @brief  Длина трубы, [м]
	double length = 80e3;
	/// @brief Высотная отметка начала трубы, [м]
	double height_in = 50;
	/// @brief Высотная отметка конца трубы, [м]
	double height_out = 100;
	/// @brief Внешний диаметр, [м]
	double external_diameter = 720e-3;
	/// @brief Толщина стенки, [м]
	double wall_thickness = 10e-3;
	/// @brief Внутренний диаметр, [м]
	double internal_diameter = external_diameter - 2 * wall_thickness;
	/// @brief Абсолютная шероховатость, [м]
	double absolute_roughness = 0.015e-3;
	/// @brief Шаг сетки
	double grid_step = length / 1000;
	/// @brief Давления в узлах сетки
	vector<double> grid_pressures = vector<double>(get_points_number(), 0);
	/// @brief Координаты узлов сетки
	vector<double> grid_coordinates = vector<double>(get_points_number(), 0);
	/// @brief Высотные отметки узлов сетки
	vector<double> grid_heights = vector<double>(get_points_number(), 0);

	/// @brief Конструктор трубы, формирующий вектор координат узлов сети
	/// @param length Длина трубы [м]
	/// @param external_diameter Внешний диаметр трубы [м]
	/// @param wall_thickness Толщина стенки, [м]
	/// @param absolute_roughness Абсолютная шероховатость, [м]
	simple_pipe_t(double length, double height_in, double height_out,
		double external_diameter, double wall_thickness, double absolute_roughness) :
		length{ length }, height_in{ height_in }, height_out{ height_out },
		external_diameter{ external_diameter }, wall_thickness{ wall_thickness }, absolute_roughness{ absolute_roughness }
	{
		recalculate_grid();
	}

	/// @brief Задать новый шаг сетки со сбросом профилей
	/// @param new_step Новый шаг
	void set_grid_step(double new_step)
	{
		grid_step = new_step;
		grid_coordinates = vector<double>(get_points_number(), 0);
		grid_pressures = vector<double>(get_points_number(), 0);
		grid_heights = vector<double>(get_points_number(), 0);
		recalculate_grid();
	}

	/// @brief Пересчитать координаты и высотные отметки (для изменения шага)
	void recalculate_grid()
	{
		size_t grid_size = get_points_number();
		for (size_t i = 0; i < grid_size; i++)
		{
			grid_coordinates[i] = i * grid_step;
			grid_heights[i] = height_in + i * (height_out - height_in) / grid_size;
		}
	}
	
	/// @brief Получить число узлов расчетной сетки
	/// @return Число узлов расчетной сетки
	size_t get_points_number()
	{
		return (size_t)(length / grid_step) + 1;
	}
};

/// @brief Условие задачи расчета изотермического течения жидкости на участке трубопровода.
/// Искомые параметры инициализированы нулями.
struct steady_flow_problem_t
{
	/// @brief Ссылка на структуру с параметрами трубы, для которой решается задача
	simple_pipe_t& pipe;
	/// @brief Давление в начале участка, [Па]
	double Pin = 5e6;
	/// @brief Высотная отметка начала участка, [м]
	double Zin = 50;
	/// @brief Давление в конце участка, [Па]
	double Pout = 0.8e6;
	/// @brief Высотная отметка конца участка, [м]
	double Zout = 100;
	/// @brief Расход жидкости, [м3/c]
	double Q = 0;
	/// @brief Плотность жидкости, [кг/м3]
	double density = 870;
	/// @brief Кинематическая вязкость, [м^2/с]
	double kinematic_viscosity = 15e-6;

	steady_flow_problem_t(simple_pipe_t& pipe, double Pin, double Pout, double Q, 
		double density, double kinematic_viscosity):
		pipe{ pipe }, Pin{ Pin }, Pout{ Pout }, Q{ Q }, density{ density }, 
		kinematic_viscosity{ kinematic_viscosity }
	{
	}

};

/// @brief Шаблон обертки для передачи функции в euler_solver
/// @tparam func_t - тип передаваемой лямбда-функции
template <typename func_t>
class function_wrapper_t
{
public:
	/// @brief Функция, для которой реализована обертка
	func_t& function;
	/// @brief Ссылка на структуру с параметрами трубы
	simple_pipe_t& pipe;
	/// @brief Ссылка на структуру с условием задачи
	steady_flow_problem_t& problem;

	/// @brief Конструктор обертки
	/// @param pipe Ссылка на структуру с параметрами трубы
	/// @param problem Ссылка на структуру с условием задачи
	/// @param function Функция, для которой реализована обертка
	function_wrapper_t(auto pipe, auto problem, auto function) :
		pipe{ pipe },
		problem{problem},
		function{function}
	{
	}

	/// @brief Получить значение функции, находящейся в обертке
	/// @param x Аргумент 
	/// @param y Аргумент
	/// @param index Индекс текущего элемента в векторе значений аргумента x
	/// @return Значение функции, находящейся в обертке
	double get_function_value(double x, double y, size_t index = 0)
	{
		return function(this->pipe, this->problem, x, y, index);
	}
};

/// @brief Шаблон для солвера, решающего ДУ 1-го порядка методом Эйлера
/// @tparam func_t 
template <typename func_t>
class euler_solver_t
{
	/// @brief Правая часть уравнения, завернутая в оболочку function_wrapper_t
	function_wrapper_t<func_t>& right_part_wrapper;

public:
	/// @brief Конструктор солвера
	/// @param right_part Правая часть уравнения, завернутая в оболочку function_wrapper_t
	euler_solver_t(function_wrapper_t<func_t>& right_part):
		right_part_wrapper{ right_part }
	{
	}

	/// @brief Решение уравнения методом Эйлера
	/// @param x Вектор аргументов. Первое значение вектора считается начальным
	/// @param y Вектор значений функции, куда будут записаны результаты. Первое значение вектора считается начальным
	/// @return Костыль для избежания ошибки компиляции
	double solve(vector<double>& x, vector<double>& y)
	{
		for (size_t i = 1; i < size(x); i++)
		{
			y[i] = y[i - 1] + (x[i] - x[i - 1]) * right_part_wrapper.get_function_value(x[i - 1], y[i - 1], i - 1);
		}
		return 0;
	}
};


/// @brief Солвер для расчета изотермического течения жидкости на участке трубопровода
class steady_flow_solver_t
{
public:
	/// @brief Поиск давления по известному давлению и расходу
	/// @param pipe Ссылка на параметры трубы
	/// @param problem Ссылка на условие задачи
	/// @param solve_Pin Флаг расчета давления в начале участка
	/// @return Рассчитанное давление, [Па]
	double solve_PQ(simple_pipe_t& pipe, steady_flow_problem_t& problem, bool solve_Pin = false)
	{
		double relative_roughness = pipe.absolute_roughness / pipe.internal_diameter;
		double speed = 4 * problem.Q / (M_PI * pow(pipe.internal_diameter, 2));
		double reynolds_number = (speed * pipe.internal_diameter) / problem.kinematic_viscosity;
		double lambda = hydraulic_resistance_isaev(reynolds_number, relative_roughness);
		if (solve_Pin)
		{
			double Pin = problem.Pout + (problem.pipe.height_in - problem.pipe.height_out) * problem.density * g +
				(lambda * pipe.length * pow(speed, 2) * problem.density) / (2 * pipe.internal_diameter);
			return Pin;
		}
		else
		{
			double Pout = problem.Pin - (problem.pipe.height_in - problem.pipe.height_out) * problem.density * g -
				lambda * (pipe.length * pow(speed, 2) * problem.density) / (2 * pipe.internal_diameter);
			return Pout;
		}

	}

	/// @brief Поиск расхода по известным давлениям методом итераций
	/// @param pipe Ссылка на параметры трубы
	/// @param problem Ссылка на условие задачи
	/// @param accuracy Точность поиска (по умолчанию 1e-3)
	/// @param iteration_max_count Предельное число итерация (по умолчанию 1000)
	/// @return Расход жидкости, [м3/с]
	double solve_PP_iteration(simple_pipe_t& pipe, steady_flow_problem_t& problem, double accuracy = 1e-3, int iteration_max_count = 1000)
	{
		double lambda_prev;
		double lambda = 0.02;
		double relative_roughness = pipe.absolute_roughness / pipe.internal_diameter;
		double speed;
		double reynolds_number;

		int iteration_count = 0;
		do
		{
			if (iteration_count > iteration_max_count)
				throw std::runtime_error("Reached maximum number of iterations");

			speed = sqrt((2 * g * pipe.internal_diameter *
				(((problem.Pin - problem.Pout) / (problem.density * g)) + problem.pipe.height_in + problem.pipe.height_out) / pipe.length) / lambda);
			reynolds_number = (speed * pipe.internal_diameter) / problem.kinematic_viscosity;
			lambda_prev = lambda;
			lambda = hydraulic_resistance_isaev(reynolds_number, relative_roughness);

			iteration_count += 1;
		} while (abs(lambda - lambda_prev) > accuracy);
		
		double Q = speed * (M_PI * pow(pipe.internal_diameter, 2)) / 4;
		return Q;
	}

	/// @brief Решение задачи PQ методом Эйлера
	/// @param pipe Ссылка на параметры трубы
	/// @param problem Ссылка на условие задачи
	vector<double> solve_PQ_euler(simple_pipe_t& pipe, steady_flow_problem_t& problem)
	{
		typedef double (*func_t)(simple_pipe_t& pipe, steady_flow_problem_t& problem, double x, double y, size_t index);
		func_t right_part
		{ (func_t)(
			[](simple_pipe_t& pipe, steady_flow_problem_t& problem, double p, double x, size_t index)
			{
				/// Лямбда-функция для расчета производной dz/dx
				auto get_height_derivative
				{ [&pipe](steady_flow_problem_t& problem, size_t index)
					{
						index = (index == 0) ? 1 : index;
						return (problem.pipe.grid_heights[index] - problem.pipe.grid_heights[index - 1])
							/ ((problem.pipe.grid_coordinates[index] - problem.pipe.grid_coordinates[index - 1]));
					} 
				};
				double relative_roughness = pipe.absolute_roughness / pipe.internal_diameter;
				double speed = 4 * problem.Q / (M_PI * pow(pipe.internal_diameter, 2));
				double reynolds_number = (speed * pipe.internal_diameter) / problem.kinematic_viscosity;
				double lambda = hydraulic_resistance_isaev(reynolds_number, relative_roughness);

				return (((-4.0) / pipe.internal_diameter) * (lambda/8.0) * problem.density * pow(speed, 2)) - 
					(problem.density * g * get_height_derivative(problem, index));
			} )
		};

		function_wrapper_t<func_t> right_part_wrapper{ pipe, problem, right_part };
		euler_solver_t<func_t> solver{ right_part_wrapper };
		pipe.grid_pressures[0] = problem.Pin;
		solver.solve(pipe.grid_coordinates, pipe.grid_pressures);
		return pipe.grid_pressures;
	}
};


/// @brief Вывод профиля давления в консоль
/// @param pipe Ссылка на параметры трубы
/// @param filename Имя файла
void pressure_to_file_(simple_pipe_t& pipe, string filename = "data.txt")
{
	std::ofstream out;
	out.open(filename);
	out << "время" << ',' << "координата" << ',' << "давление" << endl;
	if (out.is_open())
		for (size_t i = 0; i < pipe.grid_coordinates.size(); i++)
			out << 0 << ',' << pipe.grid_coordinates[i] << ',' << pipe.grid_pressures[i] << endl;
	out.close();
}


/// Оболочка для уравнения бернулли в задаче PP для использования решателя Ньютона - Рафсона
class sample_system : public fixed_system_t<1>
{
	/// @brief Ссылка на структуру с параметрами трубы
	simple_pipe_t& pipe;
	/// @brief Ссылка на структуру с условием задачи
	steady_flow_problem_t& problem;

	using fixed_system_t<1>::var_type;
public:

	/// @brief Констуктор уравнения Бернулли для задачи PP
	/// @param pipe Ссылка на сущность трубы
	/// @param problem Ссылка на сущность с условием задачи
	sample_system(simple_pipe_t& pipe, steady_flow_problem_t& problem) :
		pipe{ pipe }, problem{ problem }
	{
	}

	/// @brief Функция невязок - все члены уравнения Бернулли в правой части
	/// @param speed - скорость течения нефти, [м/с]
	/// @return Значение функци при заданной скорости
	var_type residuals(const var_type& speed) 
	{
		double relative_roughness = pipe.absolute_roughness / pipe.internal_diameter;
		double reynolds_number = (speed * pipe.internal_diameter) / problem.kinematic_viscosity;
		double lambda = hydraulic_resistance_isaev(reynolds_number, relative_roughness);
		double res = (lambda * (pipe.length * pow(speed, 2) / (pipe.internal_diameter * 2 * g))) + 
			((problem.Pout / (problem.density * g)) + problem.Zout) - 
			(((problem.Pin / (problem.density * g)) + problem.Zin));
		return { res };
	}
};

void testing_newton_raphson_()
{
	/// Задаем параметры трубы и создаем сущность трубы
	double length = 80e3;
	double external_diameter = 720e-3;
	double wall_thickness = 10e-3;
	double absolute_roughness = 0.015e-3;
	double Zin = 50;
	double Zout = 100;
	simple_pipe_t pipe{ length, Zin, Zout, external_diameter, wall_thickness, absolute_roughness };


	/// Условие для задачи PQ - расчет давления в начале участка
	double Pin = 0;
	double Pout = 0.6e6;
	double Q = 3500.0 / 3600;
	double density = 870;
	double kinematic_viscosity = 15e-6;

	/// Создание сущности, содержащей условие задачи
	steady_flow_problem_t problem1(pipe, Pin, Pout, Q, density, kinematic_viscosity);

	/// Создание солвера (не из pde_solvers)
	steady_flow_solver_t solver{};

	problem1.Pin = solver.solve_PQ(pipe, problem1, true);
	
	/// Условие для задачи PP - расчет расхода
	Pin = 5e6;
	Pout = 0.8e6;
	Q = 0;
	density = 870;
	kinematic_viscosity = 15e-6;

	/// Создание сущности, содержащей условие задачи
	steady_flow_problem_t problem3{ pipe, Pin, Pout, Q, density, kinematic_viscosity };

	/// Создание сущности решаемого уравнения
	sample_system test(pipe, problem3);

	// Задание настроек решателя по умолчанию
	fixed_solver_parameters_t<1, 0> parameters;
	// Создание структуры для записи результатов расчета
	fixed_solver_result_t<1> result;
	// Решение задачи PP с помощью решателя Ньютона - Рафсона
	fixed_newton_raphson<1>::solve_dense(test, {10}, parameters, &result);
	cout << "PP_task	Newtow_Raphson	Q = " << result.argument * (M_PI * pow(pipe.internal_diameter, 2) / 4) * 3600 << endl;
	//cout << "PP_task		Iteration	Q =  " << solver.solve_PP_iteration(pipe, problem3) * 3600 << endl;

}

void testing_PP_and_PQ_()
{
	double length = 80e3;
	double external_diameter = 720e-3;
	double wall_thickness = 10e-3;
	double absolute_roughness = 0.015e-3;
	double Zin = 50;
	double Zout = 50;
	simple_pipe_t pipe{ length, Zin, Zout, external_diameter, wall_thickness, absolute_roughness };
	steady_flow_solver_t solver{};

	/// Расчет давления в начале участка
	double Pin = 0;
	double Pout = 0.6e6;
	double Q = 3500.0/3600;
	double density = 870;
	double kinematic_viscosity = 15e-6;
	steady_flow_problem_t problem1 ( pipe, Pin, Pout, Q, density, kinematic_viscosity );
	//cout << "PoutQ_task	Formula		Pin = " << solver.solve_PQ(pipe, problem1, true) << endl;

	/// Расчет давления в конце участка
	Pin = 5.65e6;
	Pout = 0;
	Q = 3500.0 / 3600;
	density = 870;
	kinematic_viscosity = 15e-6;
	steady_flow_problem_t problem2{ pipe, Pin, Pout, Q, density, kinematic_viscosity };
	cout << "PinQ_task	Formula		Pout = " << solver.solve_PQ(pipe, problem2) << endl;
	/// Расчет задачи PQ методом Эйлера
	solver.solve_PQ_euler(pipe, problem2);
	cout << "PinQ_task	Euler		Pout = " << pipe.grid_pressures[pipe.get_points_number() - 1] << endl;
	pressure_to_file_(pipe);

	/// Расчет расхода
	Pin = 5e6;
	Pout = 0.8e6;
	Q = 0;
	density = 870;
	kinematic_viscosity = 15e-6;
	steady_flow_problem_t problem3{ pipe, Pin, Pout, Q, density, kinematic_viscosity };
	//cout << "PP_task		Iteration	Q =  " << solver.solve_PP_iteration(pipe, problem3) << endl;
}

int main()
{
	testing_newton_raphson_();
}

