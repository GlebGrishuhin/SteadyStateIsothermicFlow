#pragma once

/// @brief Формирование расчетной сетки трубы с равномерным шагом
/// @param pipe Ссылка на сущность трубы
/// @param length Длина трубы
/// @param grid_size Количество узлов сетки
void fill_profile_coordinates_(pipe_properties_t& pipe, double length, size_t grid_size)
{
    vector<double> grid(grid_size, 0);
    double grid_step = length / (grid_size - 1);
    for (size_t i = 0; i < grid_size; i++)
        grid[i] = i * grid_step;
    pipe.profile.coordinates = grid;
}

/// @brief Формирование профиля высот по высотной отметке начала и конца трубы
/// @param pipe Ссылка на сущность трубы
/// @param height_in Высотная отметка начала
/// @param height_out Высотная отметка конца
void fill_profile_heights_(pipe_properties_t& pipe, double height_in, double height_out)
{
    vector<double> heights(pipe.profile.coordinates.size(), 0);
    size_t grid_size = pipe.profile.coordinates.size();
    double height_step = (height_out - height_in) / (grid_size - 1);
    for (size_t i = 0; i < grid_size; i++)
        heights[i] = height_in + i * height_step;
    pipe.profile.heights = heights;
}

/// @brief Вывод профиля параметра в файл для построения графика скриптом на Python
/// @param pipe Ссылка на сущность трубы
/// @param profile Ссылка на профиль
/// @param filename Имя файла
void pressure_to_file_(const pipe_properties_t& pipe, const layer_t& profile, string filename = "data.txt")
{
    std::ofstream out;
    out.open(filename);
    out << "время" << ',' << "координата" << ',' << "давление" << endl;
    if (out.is_open())
        for (size_t i = 0; i < pipe.profile.coordinates.size(); i++)
            out << 0 << ',' << pipe.profile.coordinates[i] << ',' << profile[i] << endl;
    out.close();
}

void profiles_to_file_(const pipe_properties_t& pipe, const vector<vector<layer_t>>& profiles, const vector<string>& param_names, string filename = "data.txt")
{
    std::ofstream out;
    out.open(filename);
    out << "время" << ',' << "координата"; //<< ',' << "давление" << endl;
    for (size_t i = 0; i < param_names.size(); i++)
        out << ',' << param_names[i];
    out << endl;
    for (size_t time = 0; time < (profiles[0].size()); time++)
    {
        for (size_t i = 0; i < pipe.profile.coordinates.size(); i++)
        {
            out << time << ',' << pipe.profile.coordinates[i];
            for (size_t param_index = 0; param_index < param_names.size(); param_index++)
                out << ',' << profiles[param_index][time][i];
            out << endl;
        }
    }
    out.close();
}

/// @brief Задача PQ. Расчет давления в начале трубы по прямой формуле. См [Лурье 2012] раздел 4.1 задача 1
TEST(PQ_task, Formula)
{
    double length = 80e3;
    double external_diameter = 720e-3;
    double wall_thickness = 10e-3;
    double absolute_roughness = 0.015e-3;
    double Zin = 50;
    double Zout = 100;

    double density = 870;
    double kinematic_viscosity = 15e-6;

    size_t grid_size = 500;

    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    fill_profile_coordinates_(pipe, length, grid_size);
    fill_profile_heights_(pipe, Zin, Zout);

    pipe.wall.wallThickness = wall_thickness;
    pipe.wall.diameter = external_diameter - 2 * wall_thickness;
    pipe.wall.equivalent_roughness = absolute_roughness * pipe.wall.diameter;


    oil_parameters_t oil;
    oil.density.nominal_density = density;
    oil.viscosity.nominal_viscosity = kinematic_viscosity;


    /// Расчет давления в начале участка
    double Pout = 0.6e6;
    double Q = 3500.0 / 3600;
    
    double Pin = solve_PQ_(pipe, oil, Pout, Q, true);

}

/// @brief Задача PQ. Расчет профиля давлений методом Эйлера. См [Лурье 2012] раздел 4.1 задача 2
TEST(PQ_task, Euler)
{
    double length = 80e3;
    double external_diameter = 720e-3;
    double wall_thickness = 10e-3;
    double absolute_roughness = 0.015e-3;
    double Zin = 50;
    double Zout = 100;

    double density = 870;
    double kinematic_viscosity = 15e-6;


    size_t grid_size = 500;

    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    fill_profile_coordinates_(pipe, length, grid_size);
    fill_profile_heights_(pipe, Zin, Zout);

    pipe.wall.wallThickness = wall_thickness;
    pipe.wall.diameter = external_diameter - 2 * wall_thickness;
    pipe.wall.equivalent_roughness = absolute_roughness * pipe.wall.diameter;


    oil_parameters_t oil;
    oil.density.nominal_density = density;
    oil.viscosity.nominal_viscosity = kinematic_viscosity;


    /// Создаем буфер из двух слоев, каждый совпадает по размеру с pipe.profile.heights
    ring_buffer_t<vector<double>> buffer(2, vector<double>(pipe.profile.getPointCount(), 0));

    /// Задаем объемнй расход нефти, [м3/с]
    double Q = 3500.0 / 3600;

    /// Создаем расчетную модель трубы
    Pipe_model_for_PQ_t pipeModel(pipe, oil, Q);

    /// Получаем указатель на начало слоя в буфере
    profile_wrapper<double, 1> start_layer(buffer.current());

    ///Задаем начальное давление
    double Pout = 0.6e6;

    /// Модифицированный метод Эйлера для модели pipeModel,
    /// расчет ведется справа-налево относительно сетки,
    /// начальное условие Pout, 
    /// результаты расчета запишутся в слой, на который указывает start_layer
    solve_euler_corrector<1>(pipeModel, -1, Pout, &start_layer);
}

/// @brief Задача PP. Расчет расхода методом простых итераций. См [Лурье 2012] раздел 4.1 задача 2
TEST(PP_task, Iteration)
{
    double length = 80e3;
    double external_diameter = 720e-3;
    double wall_thickness = 10e-3;
    double absolute_roughness = 0.015e-3;
    double Zin = 50;
    double Zout = 100;

    double density = 870;
    double kinematic_viscosity = 15e-6;


    size_t grid_size = 500;

    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    fill_profile_coordinates_(pipe, length, grid_size);
    fill_profile_heights_(pipe, Zin, Zout);

    pipe.wall.wallThickness = wall_thickness;
    pipe.wall.diameter = external_diameter - 2 * wall_thickness;
    pipe.wall.equivalent_roughness = absolute_roughness * pipe.wall.diameter;


    oil_parameters_t oil;
    oil.density.nominal_density = density;
    oil.viscosity.nominal_viscosity = kinematic_viscosity;


    /// Расчет расхода
    double Pin = 5e6;
    double Pout = 0.8e6;

    double Q = solve_PP_iteration_(pipe, oil, Pin, Pout);

}

/// @brief Задача PP. Расчет расхода методом Ньютона-Рафсона.
TEST(PP_task, Newton)
{
    double length = 80e3;
    double external_diameter = 720e-3;
    double wall_thickness = 10e-3;
    double absolute_roughness = 0.015e-3;
    double Zin = 50;
    double Zout = 100;

    double density = 870;
    double kinematic_viscosity = 15e-6;

    size_t grid_size = 500;

    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    fill_profile_coordinates_(pipe, length, grid_size);
    fill_profile_heights_(pipe, Zin, Zout);

    pipe.wall.wallThickness = wall_thickness;
    pipe.wall.diameter = external_diameter - 2 * wall_thickness;
    pipe.wall.equivalent_roughness = absolute_roughness * pipe.wall.diameter;

    oil_parameters_t oil;
    oil.density.nominal_density = density;
    oil.viscosity.nominal_viscosity = kinematic_viscosity;

    /// Расчет расхода
    double Pin = 5e6;
    double Pout = 0.8e6;

    PP_solver_t solver(pipe, oil);
    double Q = solver.solve(Pin, Pout);
}

/// @brief Задача PP поверх Эйлера на методе Ньютона
TEST(PP_task, Euler_Newton)
{
    double length = 80e3;
    double external_diameter = 720e-3;
    double wall_thickness = 10e-3;
    double absolute_roughness = 0.015e-3;
    double Zin = 50;
    double Zout = 100;

    double density = 870;
    double kinematic_viscosity = 15e-6;

    size_t grid_size = 500;

    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    fill_profile_coordinates_(pipe, length, grid_size);
    fill_profile_heights_(pipe, Zin, Zout);

    pipe.wall.wallThickness = wall_thickness;
    pipe.wall.diameter = external_diameter - 2 * wall_thickness;
    pipe.wall.equivalent_roughness = absolute_roughness;

    oil_parameters_t oil;
    oil.density.nominal_density = density;
    oil.viscosity.nominal_viscosity = kinematic_viscosity;

    /// Расчет расхода
    double Pin = 5e6;
    double Pout = 0.8e6;

    /// Создаем буфер из двух слоев, каждый совпадает по размеру с pipe.profile.heights
    ring_buffer_t<vector<double>> buffer(2, vector<double>(pipe.profile.getPointCount(), 0));

    PP_solver_Newton_Euler_t solver(pipe, oil, buffer.current());
    double Q = solver.solve(Pin, Pout);
}