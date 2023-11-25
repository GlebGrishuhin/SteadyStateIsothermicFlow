#pragma once

typedef vector<double> layer_t;

/// @brief ������� ������ ������� ��������������� ������� �������� �� ������� ������������.
/// ������� ��������� ���������������� ������.
struct steady_flow_problem_t
{
	/// @brief ������ �� ��������� � ����������� �����, ��� ������� �������� ������
	const pipe_properties_t& pipe;
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

	steady_flow_problem_t(pipe_properties_t& pipe, double Pin, double Pout, double Q,
		double density, double kinematic_viscosity) :
		pipe{ pipe }, Pin{ Pin }, Pout{ Pout }, Q{ Q }, density{ density },
		kinematic_viscosity{ kinematic_viscosity }
	{
	}

};

/// @brief ������ ������� ��� �������� ������� � euler_solver
/// @tparam func_t - ��� ������������ ������-�������
template <typename func_t>
class function_wrapper_t
{
public:
	/// @brief �������, ��� ������� ����������� �������
	func_t& function;
	/// @brief ������ �� ��������� � ����������� �����
	pipe_properties_t& pipe;
	/// @brief ������ �� ��������� � �������� ������
	steady_flow_problem_t& problem;

	/// @brief ����������� �������
	/// @param pipe ������ �� ��������� � ����������� �����
	/// @param problem ������ �� ��������� � �������� ������
	/// @param function �������, ��� ������� ����������� �������
	function_wrapper_t(auto pipe, auto problem, auto function) :
		pipe{ pipe },
		problem{ problem },
		function{ function }
	{
	}

	/// @brief �������� �������� �������, ����������� � �������
	/// @param x �������� 
	/// @param y ��������
	/// @param index ������ �������� �������� � ������� �������� ��������� x
	/// @return �������� �������, ����������� � �������
	double get_function_value(double x, double y, size_t index = 0)
	{
		return function(this->pipe, this->problem, x, y, index);
	}
};

/// @brief ������ ��� �������, ��������� �� 1-�� ������� ������� ������
/// @tparam func_t 
template <typename func_t>
class euler_solver_t
{
	/// @brief ������ ����� ���������, ���������� � �������� function_wrapper_t
	function_wrapper_t<func_t>& right_part_wrapper;

public:
	/// @brief ����������� �������
	/// @param right_part ������ ����� ���������, ���������� � �������� function_wrapper_t
	euler_solver_t(function_wrapper_t<func_t>& right_part) :
		right_part_wrapper{ right_part }
	{
	}

	/// @brief ������� ��������� ������� ������
	/// @param x ������ ����������. ������ �������� ������� ��������� ���������
	/// @param y ������ �������� �������, ���� ����� �������� ����������. ������ �������� ������� ��������� ���������
	/// @return ������� ��� ��������� ������ ����������
	double solve(vector<double>& x, vector<double>& y)
	{
		for (size_t i = 1; i < size(x); i++)
		{
			y[i] = y[i - 1] + (x[i] - x[i - 1]) * right_part_wrapper.get_function_value(x[i - 1], y[i - 1], i - 1);
		}
		return 0;
	}
};


/// @brief ����� �������� �� ���������� �������� � �������
/// @param pipe ������ �� ��������� �����
/// @param problem ������ �� ������� ������
/// @param solve_Pin ���� ������� �������� � ������ �������
/// @return ������������ ��������, [��]
double solve_PQ_(const pipe_properties_t& pipe, const oil_parameters_t& oil, double* init_cond, bool solve_Pin = false)
{
	// ��������� ��������� ������ ������
	double P = init_cond[0];
	double Q = init_cond[1];

	// ��������� �����
	// ���������� ������� ��� �������� � ������� ������ �������
	double rho = oil.density();
	double mu = oil.viscosity();
	// ��������� �����
	double relative_roughness = pipe.wall.relativeRoughness();
	double d = pipe.wall.diameter;
	double z0 = pipe.profile.heights.front();
	double zL = pipe.profile.heights.back();
	double L = pipe.profile.getLength();

	double v = 4 * Q / (M_PI * pow(d, 2));
	double reynolds_number = (v * d) / mu;;
	double lambda = pipe.resistance_function(reynolds_number, relative_roughness);

	// ������ ������ (������ ����� ��������� ��������)
	double dH = lambda * L * pow(v, 2) / (2 * d * M_G);
	// ������� �����
	double dZ = zL - z0;
	if (solve_Pin)
	{
		double HL = P / (rho * M_G);
		double Pin = (HL + dH + dZ) * rho * M_G;
		return Pin;
	}
	else
	{
		double H0 = P / (rho * M_G);
		double Pout = (H0 - dH - dZ) * rho * M_G;
		return Pout;
	}
}



/// @brief ����� ������� �� ��������� ��������� ������� ��������
/// @param pipe ������ �� ��������� �����
/// @param problem ������ �� ������� ������
/// @param accuracy �������� ������ (�� ��������� 1e-3)
/// @param iteration_max_count ���������� ����� �������� (�� ��������� 1000)
/// @return ������ ��������, [�3/�]
double solve_PP_iteration_(const pipe_properties_t& pipe, const oil_parameters_t& oil, double pressures[2],
	double accuracy = 1e-4, int iteration_max_count = 1000)
{
	double lambda_prev;
	// ��������� ����������� ������
	double lambda = 0.02;
	// �������� ������
	double v;
	// ����� ����������
	double reynolds_number;
	// ��������� ������� ������ �� �����
	double dH;

	// ��������� �����
	// ��������� ���������� ������� ��� �������� � ������� ������ ������ �������
	double rho = oil.density();
	double mu = oil.viscosity();
	// ��������� �����
	double z0 = pipe.profile.heights.front();
	double zL = pipe.profile.heights.back();
	double L = pipe.profile.getLength();
	double d = pipe.wall.diameter;
	double relative_roughness = pipe.wall.relativeRoughness();

	double Pin = pressures[0];
	double Pout = pressures[1];

	int iteration_count = 0;
	do
	{
		if (iteration_count > iteration_max_count)
			throw std::runtime_error("Reached maximum number of iterations");

		dH = (Pin - Pout) / (rho * M_G) + (z0 - zL);
		v = sqrt(2 * M_G * d * dH / (L * lambda));
		reynolds_number = (v * d) / mu;
		lambda_prev = lambda;
		lambda = pipe.resistance_function(reynolds_number, relative_roughness);
		iteration_count += 1;
	} while (abs(lambda - lambda_prev) > accuracy);
	double Q = v * (M_PI * pow(d, 2)) / 4;
	return Q;
}


/// @brief ����� ������� �������� � �������
/// @param pipe ������ �� ��������� �����
/// @param filename ��� �����
void pressure_to_file_(const pipe_properties_t& pipe, const layer_t& profile, string filename = "data.txt")
{
	std::ofstream out;
	out.open(filename);
	out << "�����" << ',' << "����������" << ',' << "��������" << endl;
	if (out.is_open())
		for (size_t i = 0; i < pipe.profile.coordinates.size(); i++)
			out << 0 << ',' << pipe.profile.coordinates[i] << ',' << profile[i] << endl;
	out.close();
}


/// �������� ��� ��������� �������� � ������ PP ��� ������������� �������� ������� - �������
class PP_solver_t : public fixed_system_t<1>
{
	/// @brief ������ �� ��������� � ����������� �����
	const oil_parameters_t& oil;
	/// @brief ������ �� ��������� � �������� ������
	const pipe_properties_t& pipe;
	/// @brief ������� ��������, �� �������� ����� ����� �������� � ������ � ����� �����
	layer_t& pressure_profile;

	using fixed_system_t<1>::var_type;
public:

	/// @brief ���������� ��������� �������� ��� ������ PP
	/// @param pipe ������ �� �������� �����
	/// @param problem ������ �� �������� � �������� ������
	PP_solver_t(const pipe_properties_t& pipe, const oil_parameters_t& oil, layer_t pressure_profile) :
		pipe{ pipe }, oil{ oil }, pressure_profile{ pressure_profile }
	{
	}

	/// @brief ������� ������� - ��� ����� ��������� �������� � ������ �����
	/// @param speed - �������� ������� �����, [�/�]
	/// @return �������� ������� ������� ��� �������� ��������
	var_type residuals(const var_type& v)
	{
		// ��������� ����� � �����
		// ��������� ���������� ������� ��� �������� � ������� ������ ������ �������
		double rho = oil.density();
		double mu = oil.viscosity();

		double z0 = pipe.profile.heights.front();
		double zL = pipe.profile.heights.back();
		double L = pipe.profile.getLength();
		double d = pipe.wall.diameter;
		double relative_roughness = pipe.wall.relativeRoughness();

		double p0 = pressure_profile.front();
		double pL = pressure_profile.back();

		double reynolds_number = (v * d) / mu;
		double lambda = pipe.resistance_function(reynolds_number, relative_roughness);

		double H0 = p0 / (rho * M_G) + z0;
		double HL = pL / (rho * M_G) + zL;
		double dH = lambda * L * pow(v, 2) / (d * 2 * M_G);
		double result = H0 - HL - dH;
		return result;
	}

	double solve(const double initial_argument = 0.1)
	{
		// ������� �������� �������� �� ���������
		fixed_solver_parameters_t<1, 0> parameters;
		// �������� ��������� ��� ������ ����������� �������
		fixed_solver_result_t<1> result;
		// ������� ������ PP � ������� �������� ������� - �������
		fixed_newton_raphson<1>::solve_dense(*this, initial_argument, parameters, &result);
		return result.argument;
	}

};

void testing_newton_raphson_()
{

	oil_parameters_t oil;
	pipe_properties_t pipe;


	/// ������ ��������� ����� � ������� �������� �����
	double length = 80e3;
	double external_diameter = 720e-3;
	double wall_thickness = 10e-3;
	double absolute_roughness = 0.015e-3;
	double Zin = 50;
	double Zout = 100;
	pipe_properties_t pipe{ length, Zin, Zout, external_diameter, wall_thickness, absolute_roughness };


	/// ������� ��� ������ PQ - ������ �������� � ������ �������
	double Pin = 0;
	double Pout = 0.6e6;
	double Q = 3500.0 / 3600;
	double density = 870;
	double kinematic_viscosity = 15e-6;

	/// �������� ��������, ���������� ������� ������
	steady_flow_problem_t problem1(pipe, Pin, Pout, Q, density, kinematic_viscosity);

	/// �������� ������� (�� �� pde_solvers)
	steady_flow_solver_t solver{};

	problem1.Pin = solver.solve_PQ(pipe, problem1, true);

	/// ������� ��� ������ PP - ������ �������
	Pin = 5e6;
	Pout = 0.8e6;
	Q = 0;
	density = 870;
	kinematic_viscosity = 15e-6;

	/// �������� ��������, ���������� ������� ������
	steady_flow_problem_t problem3{ pipe, Pin, Pout, Q, density, kinematic_viscosity };

	/// �������� �������� ��������� ���������
	sample_system test(pipe, problem3);

	// ������� �������� �������� �� ���������
	fixed_solver_parameters_t<1, 0> parameters;
	// �������� ��������� ��� ������ ����������� �������
	fixed_solver_result_t<1> result;
	// ������� ������ PP � ������� �������� ������� - �������
	fixed_newton_raphson<1>::solve_dense(test, 10, parameters, &result);
	cout << "PP_task	Newtow_Raphson	Q = " << result.argument * (M_PI * pow(pipe.internal_diameter, 2) / 4) * 3600 << " m3/h" << endl;
	//cout << "PP_task		Iteration	Q =  " << solver.solve_PP_iteration(pipe, problem3) * 3600 << " m3/h" << endl;

}

void testing_PP_and_PQ_()
{
	double length = 80e3;
	double external_diameter = 720e-3;
	double wall_thickness = 10e-3;
	double absolute_roughness = 0.015e-3;
	double Zin = 50;
	double Zout = 100;
	pipe_properties_t pipe{ length, Zin, Zout, external_diameter, wall_thickness, absolute_roughness };
	steady_flow_solver_t solver{};

	/// ������ �������� � ������ �������
	double Pin = 0;
	double Pout = 0.6e6;
	double Q = 3500.0 / 3600;
	double density = 870;
	double kinematic_viscosity = 15e-6;
	steady_flow_problem_t problem1(pipe, Pin, Pout, Q, density, kinematic_viscosity);
	//cout << "PoutQ_task	Formula		Pin = " << solver.solve_PQ(pipe, problem1, true) << endl;

	/// ������ �������� � ����� �������
	Pin = 5.65e6;
	Pout = 0;
	Q = 3500.0 / 3600;
	density = 870;
	kinematic_viscosity = 15e-6;
	steady_flow_problem_t problem2{ pipe, Pin, Pout, Q, density, kinematic_viscosity };
	cout << "PinQ_task	Formula		Pout = " << solver.solve_PQ(pipe, problem2) << endl;
	/// ������ ������ PQ ������� ������
	solver.solve_PQ_euler(pipe, problem2);
	cout << "PinQ_task	Euler		Pout = " << pipe.grid_pressures[pipe.get_points_number() - 1] << endl;
	pressure_to_file_(pipe);

	/// ������ �������
	Pin = 5e6;
	Pout = 0.8e6;
	Q = 0;
	density = 870;
	kinematic_viscosity = 15e-6;
	steady_flow_problem_t problem3{ pipe, Pin, Pout, Q, density, kinematic_viscosity };
	cout << "PP_task		Iteration	Q =  " << solver.solve_PP_iteration(pipe, problem3) << endl;
}

int main()
{
	//testing_newton_raphson_();
	testing_PP_and_PQ_();
}

