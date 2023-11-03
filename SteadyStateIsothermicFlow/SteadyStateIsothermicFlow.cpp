#include <iomanip>
#include <iostream>


#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

using namespace std;

const double g = 9.81;

/// @brief Физические свойства трубы
struct simple_pipe_t
{
	/// @brief  Длина трубы, [м]
	double length = 80e3;
	/// @brief Внешний диаметр, [м]
	double external_diameter = 720e-3;
	/// @brief Толщина стенки, [м]
	double wall_thickness = 10e-3;
	/// @brief Внутренний диаметр, [м]
	double internal_diameter = external_diameter - 2 * wall_thickness;
	/// @brief Абсолютная шероховатость, [м]
	double absolute_roughness = 0.015e-3;
};

/// @brief Условие задачи расчета изотермического течения жидкости на участке трубопровода.
/// Искомые параметры инициализированы нулями.
struct steady_flow_problem_t
{
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
			double Pin = problem.Pout + (problem.Zin - problem.Zout) * problem.density * g + 
				(lambda * pipe.length * pow(speed, 2) * problem.density) / (2 * pipe.internal_diameter);
			return Pin;
		}
		else
		{
			double Pout = problem.Pin - (problem.Zin - problem.Zout) * problem.density * g - 
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
		double lambda_prev = -INFINITY;
		double lambda = 0.02;
		double relative_roughness = pipe.absolute_roughness / pipe.internal_diameter;
		double speed;
		double reynolds_number;

		int iteration_count = 0;

		while (abs(lambda - lambda_prev) > accuracy)
		{
			if (iteration_count > iteration_max_count)
				throw std::runtime_error("Reached maximum number of iterations");

			speed = sqrt((2 * g * pipe.internal_diameter * 
				(((problem.Pin - problem.Pout) / (problem.density * g)) + problem.Zin - problem.Zout) / pipe.length) / lambda);
			reynolds_number = (speed * pipe.internal_diameter) / problem.kinematic_viscosity;
			lambda_prev = lambda;
			lambda = hydraulic_resistance_isaev(reynolds_number, relative_roughness);
			
			iteration_count += 1;
		}
		double Q = speed * (M_PI * pow(pipe.internal_diameter, 2)) / 4;
		return Q;
	}
};

int main()
{
	simple_pipe_t pipe{};
	steady_flow_solver_t solver{};

	/// Расчет давления в начале участка
	double Pin = 0;
	double Zin = 50;
	double Pout = 0.6e6;
	double Zout = 100;
	double Q = 3500.0/3600;
	double density = 870;
	double kinematic_viscosity = 15e-6;
	steady_flow_problem_t problem{ Pin, Zin, Pout, Zout, Q, density, kinematic_viscosity };
	cout << "Pin = " << solver.solve_PQ(pipe, problem, true) << endl;

	/// Расчет давления в конце участка
	Pin = 5.65e6;
	Zin = 50;
	Pout = 0;
	Zout = 100;
	Q = 3500.0 / 3600;
	density = 870;
	kinematic_viscosity = 15e-6;
	problem = steady_flow_problem_t{ Pin, Zin, Pout, Zout, Q, density, kinematic_viscosity };
	cout << "Pout " << solver.solve_PQ(pipe, problem) << endl;

	/// Расчет расхода
	Pin = 5e6;
	Zin = 50;
	Pout = 0.8e6;
	Zout = 100;
	Q = 0;
	density = 870;
	kinematic_viscosity = 15e-6;
	problem = steady_flow_problem_t{ Pin, Zin, Pout, Zout, Q, density, kinematic_viscosity };
	cout << "Q =  " << solver.solve_PP_iteration(pipe, problem) << endl;
}

