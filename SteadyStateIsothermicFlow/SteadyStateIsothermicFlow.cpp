#include <iomanip>
#include <iostream>


#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

using namespace std;

const double g = 9.81;

/// @brief ���������� �������� �����
struct simple_pipe_t
{
	/// @brief  ����� �����, [�]
	double length = 80e3;
	/// @brief ������� �������, [�]
	double external_diameter = 720e-3;
	/// @brief ������� ������, [�]
	double wall_thickness = 10e-3;
	/// @brief ���������� �������, [�]
	double internal_diameter = external_diameter - 2 * wall_thickness;
	/// @brief ���������� �������������, [�]
	double absolute_roughness = 0.015e-3;
};

/// @brief ������� ������ ������� ��������������� ������� �������� �� ������� ������������.
/// ������� ��������� ���������������� ������.
struct steady_flow_problem_t
{
	/// @brief �������� � ������ �������, [��]
	double Pin = 5e6;
	/// @brief �������� ������� ������ �������, [�]
	double Zin = 50;
	/// @brief �������� � ����� �������, [��]
	double Pout = 0.8e6;
	/// @brief �������� ������� ����� �������, [�]
	double Zout = 100;
	/// @brief ������ ��������, [�3/c]
	double Q = 0;
	/// @brief ��������� ��������, [��/�3]
	double density = 870;
	/// @brief �������������� ��������, [�^2/�]
	double kinematic_viscosity = 15e-6;
};

/// @brief ������ ��� ������� ��������������� ������� �������� �� ������� ������������
class steady_flow_solver_t
{
public:
	/// @brief ����� �������� �� ���������� �������� � �������
	/// @param pipe ������ �� ��������� �����
	/// @param problem ������ �� ������� ������
	/// @param solve_Pin ���� ������� �������� � ������ �������
	/// @return ������������ ��������, [��]
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

	/// @brief ����� ������� �� ��������� ��������� ������� ��������
	/// @param pipe ������ �� ��������� �����
	/// @param problem ������ �� ������� ������
	/// @param accuracy �������� ������ (�� ��������� 1e-3)
	/// @param iteration_max_count ���������� ����� �������� (�� ��������� 1000)
	/// @return ������ ��������, [�3/�]
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

	/// ������ �������� � ������ �������
	double Pin = 0;
	double Zin = 50;
	double Pout = 0.6e6;
	double Zout = 100;
	double Q = 3500.0/3600;
	double density = 870;
	double kinematic_viscosity = 15e-6;
	steady_flow_problem_t problem{ Pin, Zin, Pout, Zout, Q, density, kinematic_viscosity };
	cout << "Pin = " << solver.solve_PQ(pipe, problem, true) << endl;

	/// ������ �������� � ����� �������
	Pin = 5.65e6;
	Zin = 50;
	Pout = 0;
	Zout = 100;
	Q = 3500.0 / 3600;
	density = 870;
	kinematic_viscosity = 15e-6;
	problem = steady_flow_problem_t{ Pin, Zin, Pout, Zout, Q, density, kinematic_viscosity };
	cout << "Pout " << solver.solve_PQ(pipe, problem) << endl;

	/// ������ �������
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

