﻿#pragma once

#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>


using namespace std;

using namespace pde_solvers;

typedef vector<double> layer_t;

/// @brief Поиск давления по известному давлению и расходу. См[Лурье 2012] раздел 4.1 задача 1
/// @param pipe Ссылка на сущность трубы
/// @param oil Ссылка на сущность нефти
/// @param P Давление в трубе, [Па]
/// @param Q Объемный расход нефти, [м3/с]
/// @param solve_Pin Флаг расчета давления в начале участка (при 1 значение P трактуется как давление в начале)
/// @return Рассчитанное давление, [Па]
double solve_PQ_(const pipe_properties_t& pipe, const oil_parameters_t& oil, double P, double Q, bool solve_Pin = false)
{
	// Параметры нефти
	// Переменные введены для удобства и краткой записи расчета
	double rho = oil.density();
	double mu = oil.viscosity();
	// Параметры трубы
	double relative_roughness = pipe.wall.relativeRoughness();
	double d = pipe.wall.diameter;
	double z0 = pipe.profile.heights.front();
	double zL = pipe.profile.heights.back();
	double L = pipe.profile.getLength();

	double v = 4 * Q / (M_PI * pow(d, 2));
	double reynolds_number = (v * d) / mu;;
	double lambda = pipe.resistance_function(reynolds_number, relative_roughness);

	// Потери напора (правая часть уравнения Бернулли)
	double dH = lambda * L * pow(v, 2) / (2 * d * M_G);
	// Перепад высот
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

/// @brief Поиск расхода по известным давлениям методом итераций
/// @param pipe Ссылка на сущность трубы
/// @param oil Ссылка на сущность нефти
/// @param Pin Давление в начале трубы, [Па]
/// @param Pout Давление в конце трубы, [Па]
/// @param accuracy Точность поиска (по умолчанию 1e-3)
/// @param iteration_max_count Предельное число итерация (по умолчанию 1000)
/// @return Расход жидкости, [м3/с]
double solve_PP_iteration_(const pipe_properties_t& pipe, const oil_parameters_t& oil, 
	double Pin, double Pout, double accuracy = 1e-4, int iteration_max_count = 1000)
{
	double lambda_prev;
	// Начальное приближение лямбды
	double lambda = 0.02;
	// Скорость потока
	double v;
	// Число Рейнольдса
	double reynolds_number;
	// Изменение полного напора по трубе
	double dH;

	// Параметры нефти
	// Отдельные переменные введены для удобства и краткой записи формул расчета
	double rho = oil.density();
	double mu = oil.viscosity();
	// Параметры трубы
	double z0 = pipe.profile.heights.front();
	double zL = pipe.profile.heights.back();
	double L = pipe.profile.getLength();
	double d = pipe.wall.diameter;
	double relative_roughness = pipe.wall.relativeRoughness();

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

/// @brief Уравнение трубы для решения задачи PQ методом Эйлера  
class Pipe_model_for_PQ_t : public ode_t<1>
{
public:
	using ode_t<1>::equation_coeffs_type;
	using ode_t<1>::right_party_type;
	using ode_t<1>::var_type;
protected:
	const pipe_properties_t& pipe;
	const oil_parameters_t& oil;
	const double flow;

public:
	/// @brief Констуктор уравнения трубы
	/// @param pipe Ссылка на сущность трубы
	/// @param oil Ссылка на сущность нефти
	/// @param flow Объемный расход, [м3/с]
	Pipe_model_for_PQ_t(const pipe_properties_t& pipe, const oil_parameters_t& oil, double flow)
		: pipe(pipe)
		, oil(oil)
		, flow(flow)
	{
	}

	/// @brief Возвращает известную уравнению сетку
	virtual const vector<double>& get_grid() const override {
		return pipe.profile.coordinates;
	}

	/// @brief Возвращает значение правой части ДУ
	/// см. файл 2023-11-09 Реализация стационарных моделей с прицелом на квазистационар.docx
	/// @param grid_index Обсчитываемый индекс расчетной сетки
	/// @param point_vector Начальные условия
	/// @return Значение правой части ДУ в точке point_vector
	virtual right_party_type ode_right_party(
		size_t grid_index, const var_type& point_vector) const override
	{
		double rho = oil.density();
		double S_0 = pipe.wall.getArea();
		double v = flow / (S_0);
		double Re = v * pipe.wall.diameter / oil.viscosity();
		double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());
		double tau_w = lambda / 8 * rho * v * abs(v);

		// Обработка индекса в случае расчетов на границах трубы
		// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
		grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
		grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

		// Расчет производной высотного профиля по координате dz/dx
		double height_derivative = (pipe.profile.heights[grid_index] - pipe.profile.heights[grid_index - 1]) /
			(pipe.profile.coordinates[grid_index] - pipe.profile.coordinates[grid_index - 1]);

		return { ((-4) / pipe.wall.diameter) * tau_w - rho * M_G * height_derivative };
	}

};

/// @brief Оболочка для уравнения Бернулли в задаче PP для использования решателя Ньютона - Рафсона
class PP_solver_t : public fixed_system_t<1>
{
	/// @brief Ссылка на структуру с параметрами трубы
	const oil_parameters_t& oil;
	/// @brief Ссылка на структуру с условием задачи
	const pipe_properties_t& pipe;
	/// @brief Поле класса для хранения давления в начале трубы
	double p0 = 0;
	/// @brief Поле класса для хранения давления в конце трубы
	double pL = 0;

	using fixed_system_t<1>::var_type;
public:

	/// @brief Констуктор уравнения Бернулли для задачи PP
	/// @param pipe Ссылка на сущность трубы
	/// @param problem Ссылка на сущность с условием задачи
	PP_solver_t(const pipe_properties_t& pipe, const oil_parameters_t& oil) :
		pipe{ pipe }, oil{ oil }
	{
	}

	/// @brief Функция невязок - все члены уравнения Бернулли в правой части
	/// @param v - скорость течения нефти, [м/с]
	/// @return Значение функции невязок при заданной скорости
	var_type residuals(const var_type& v)
	{
		// Параметры трубы и нефти
		// Отдельные переменные введены для удобства и краткой записи формул расчета
		double rho = oil.density();
		double mu = oil.viscosity();

		double z0 = pipe.profile.heights.front();
		double zL = pipe.profile.heights.back();
		double L = pipe.profile.getLength();
		double d = pipe.wall.diameter;
		double relative_roughness = pipe.wall.relativeRoughness();

		double p0 = this->p0;
		double pL = this->pL;

		double reynolds_number = (v * d) / mu;
		double lambda = pipe.resistance_function(reynolds_number, relative_roughness);

		double H0 = p0 / (rho * M_G) + z0;
		double HL = pL / (rho * M_G) + zL;
		double dH = lambda * L * pow(v, 2) / (d * 2 * M_G);
		double result = H0 - HL - dH;
		return result;
	}

	/// @brief Поиск скорости течения 
	/// @param initial_argument - начальное приближение, [м/с]
	/// @return Скорость течения, [м/с]
	double solve(double Pin, double Pout, const double initial_argument = 0.1)
	{
		this->p0 = Pin;
		this->pL = Pout;
		// Задание настроек решателя по умолчанию
		fixed_solver_parameters_t<1, 0> parameters;
		// Создание структуры для записи результатов расчета
		fixed_solver_result_t<1> result;
		// Решение задачи PP с помощью решателя Ньютона - Рафсона
		fixed_newton_raphson<1>::solve_dense(*this, initial_argument, parameters, &result);
		return result.argument * M_PI * pow(pipe.wall.diameter, 2) / 4;
	}
};


/// @brief Оболочка для уравнения Бернулли в задаче PP поверх Эйлера на методе Ньютона
class PP_solver_Newton_Euler_t : public fixed_system_t<1>
{
	/// @brief Ссылка на структуру с параметрами трубы
	const oil_parameters_t& oil;
	/// @brief Ссылка на структуру с условием задачи
	const pipe_properties_t& pipe;
	/// @brief Эталонное значение давления на входе
	double Pin_ideal = 0;
	/// @brief Эталонное значение давления на выходе
	double Pout_ideal = 0;
	/// @brief Слой буфера в буфере для записи результата расчета
	layer_t layer;

	using fixed_system_t<1>::var_type;
public:

	/// @brief Констуктор уравнения Бернулли для задачи PP
	/// @param pipe Ссылка на сущность трубы
	/// @param problem Ссылка на сущность с условием задачи
	PP_solver_Newton_Euler_t(const pipe_properties_t& pipe, const oil_parameters_t& oil, layer_t layer) :
		pipe{ pipe }, oil{ oil }, layer{ layer}
	{
	}

	/// @brief Функция невязок - отличие P0 от P0_ideal при заданной скорости потока v
	/// @param v - скорость течения нефти, [м/с]
	/// @return Значение функции невязок при заданной скорости
	var_type residuals(const var_type& v)
	{
		/// Объемный расход нефти на основе скорости потока, [м3/с]
		double Q = v * M_PI * pow(pipe.wall.diameter, 2) / 4;

		/// Создаем расчетную модель трубы
		Pipe_model_for_PQ_t pipeModel(this->pipe, this->oil, Q);

		/// Получаем указатель на начало слоя в буфере
		profile_wrapper<double, 1> start_layer(layer);

		/// Модифицированный метод Эйлера для модели pipeModel,
		/// расчет ведется справа-налево относительно сетки,
		/// начальное условие Pout_ideal, 
		/// результаты расчета запишутся в слой, на который указывает start_layer
		solve_euler_corrector<1>(pipeModel, -1, Pout_ideal, &start_layer);

		double result = layer.front() - this->Pin_ideal;
		return result;
	}

	/// @brief Поиск скорости течения 
	/// @param Pin_ideal Эталонное давление в начале трубы [Па]
	/// @param Pout_ideal Эталонное давление в конце трубы [Па]
	/// @param initial_argument - начальное приближение, [м/с]
	/// @return Скорость течения, [м/с]
	double solve(const double Pin_ideal, double Pout_ideal, const double initial_speed = 0.1)
	{
		this->Pin_ideal = Pin_ideal;
		this->Pout_ideal = Pout_ideal;
		// Задание настроек решателя по умолчанию
		fixed_solver_parameters_t<1, 0> parameters;
		// Создание структуры для записи результатов расчета
		fixed_solver_result_t<1> result;
		// Решение задачи PP с помощью решателя Ньютона - Рафсона
		fixed_newton_raphson<1>::solve_dense(*this, initial_speed, parameters, &result);
		return result.argument * M_PI * pow(pipe.wall.diameter, 2) / 4;
	}
};


/// @brief Изотермический квазистационарный гидравлический расчет для партий
class Pipe_model_oil_parties_t : public ode_t<1>
{
public:
	using ode_t<1>::equation_coeffs_type;
	using ode_t<1>::right_party_type;
	using ode_t<1>::var_type;
protected:
	const pipe_properties_t& pipe;
	const double flow;
	const layer_t& density_layer;
	const layer_t& viscosity_layer;

public:
	/// @brief Констуктор уравнения трубы
	/// @param pipe Ссылка на сущность трубы
	/// @param density_layer Ссылка на профиль плотности
	/// @param viscosity_layer Ссылка на профиль вязкости
	/// @param flow Объемный расход, [м3/с]
	Pipe_model_oil_parties_t(const pipe_properties_t& pipe,
		const layer_t& density_layer, const layer_t& viscosity_layer, double flow)
		: pipe(pipe)
		, flow(flow)
		, density_layer{ density_layer }
		, viscosity_layer{ viscosity_layer }
	{
	}

	/// @brief Возвращает известную уравнению сетку
	virtual const vector<double>& get_grid() const override {
		return pipe.profile.coordinates;
	}

	/// @brief Возвращает значение правой части ДУ
	/// см. файл 2023-11-09 Реализация стационарных моделей с прицелом на квазистационар.docx
	/// @param grid_index Обсчитываемый индекс расчетной сетки
	/// @param point_vector Начальные условия
	/// @return Значение правой части ДУ в точке point_vector
	virtual right_party_type ode_right_party(
		size_t grid_index, const var_type& point_vector) const override
	{
		
		double pressure = point_vector;
		double density = density_layer[grid_index];
		double S_0 = pipe.wall.getArea();
		double v = flow / (S_0);
		double viscosity = viscosity_layer[grid_index];
		double Re = v * pipe.wall.diameter / viscosity;
		double lambda = pipe.resistance_function(Re, pipe.wall.relativeRoughness());

		// Обработка индекса в случае расчетов на границах трубы
		// Чтобы не выйти за массив высот, будем считать dz/dx в соседней точке
		grid_index = grid_index == 0 ? grid_index + 1 : grid_index;
		grid_index = grid_index == pipe.profile.heights.size() - 1 ? grid_index - 1 : grid_index;

		// Расчет производной высотного профиля по координате dz/dx
		double height_derivative = (pipe.profile.heights[grid_index] - pipe.profile.heights[grid_index - 1]) /
			(pipe.profile.coordinates[grid_index] - pipe.profile.coordinates[grid_index - 1]);

		double density_derivative = (density_layer[grid_index] - density_layer[grid_index - 1]) /
			(pipe.profile.coordinates[grid_index] - pipe.profile.coordinates[grid_index - 1]);

		return (-1 * lambda * (1.0 / pipe.wall.diameter) * density * pow(v, 2) / 2) + height_derivative * density * M_G;
	}

};