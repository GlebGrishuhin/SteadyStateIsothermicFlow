#pragma once

void fill_profile_coordinates_(pipe_properties_t& pipe, double length, size_t grid_size)
{
    vector<double> grid(grid_size, 0);
    double grid_step = length / (grid_size - 1);
    for (size_t i = 0; i < grid_size; i++)
        grid[i] = i * grid_step;
    pipe.profile.coordinates = grid;

}

void fill_profile_heights_(pipe_properties_t& pipe, double height_in, double height_out)
{
    vector<double> heights(pipe.profile.coordinates.size(), 0);
    size_t grid_size = pipe.profile.coordinates.size();
    double height_step = (height_out - height_in) / (grid_size - 1);
    for (size_t i = 0; i < grid_size; i++)
        heights[i] = height_in + i * height_step;
    pipe.profile.heights = heights;

}

TEST(PQ_task, TemplatedLayer1)
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
    double Pin = 0;
    double Pout = 0.6e6;
    double Q = 3500.0 / 3600;
    
    solve_PQ_(pipe, oil, { Pout, Q }, true);

}