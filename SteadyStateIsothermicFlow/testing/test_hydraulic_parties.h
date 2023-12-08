#pragma once

/// @brief Решение уравнения Бернулли с учетом партий, раcсчитанных методом характеристик
TEST(Hydraulic_Parties, MOC_Euler)
{
    double length = 80e3;
    double external_diameter = 720e-3;
    double wall_thickness = 10e-3;
    double absolute_roughness = 0.015e-3;
    double Zin = 50;
    double Zout = 100;

    double Pin = 5e6;

    size_t grid_size = 500;

    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    fill_profile_coordinates_(pipe, length, grid_size);
    fill_profile_heights_(pipe, Zin, Zout);

    pipe.wall.wallThickness = wall_thickness;
    pipe.wall.diameter = external_diameter - 2 * wall_thickness;
    pipe.wall.equivalent_roughness = absolute_roughness;

    typedef composite_layer_t<profile_collection_t<1>,
        moc_solver<1>::specific_layer> single_var_moc_t;


    // Буферы для плотности, вязкости, давления
    ring_buffer_t<single_var_moc_t> buffer_density(2, pipe.profile.getPointCount());
    ring_buffer_t<single_var_moc_t> buffer_viscosity(2, pipe.profile.getPointCount());
    ring_buffer_t<single_var_moc_t> buffer_pressure(2, pipe.profile.getPointCount());

    buffer_density.advance(+1);
    buffer_viscosity.advance(+1);
    buffer_pressure.advance(+1);

    // Ссылки на текущий и предыдущий слои в буфера
    single_var_moc_t& density_layer_prev = buffer_density.previous();
    single_var_moc_t& density_layer_next = buffer_density.current();

    single_var_moc_t& viscosity_layer_prev = buffer_viscosity.previous();
    single_var_moc_t& viscosity_layer_next = buffer_viscosity.current();

    single_var_moc_t& pressure_layer_prev = buffer_pressure.previous();
    single_var_moc_t& pressure_layer_next = buffer_pressure.current();


    auto& rho_initial = density_layer_prev.vars.point_double[0];
    rho_initial = vector<double>(rho_initial.size(), 850); // инициализация начальной плотности

    auto& viscosity_initial = viscosity_layer_prev.vars.point_double[0];
    viscosity_initial = vector<double>(viscosity_initial.size(), 10e-6); // инициализация начальной вязкости

    vector<double> Q(pipe.profile.getPointCount(), 0.8); // задаем по трубе расход
    PipeQAdvection advection_model(pipe, Q);

    double rho_in = 860; // плотность нефти, закачиваемой на входе трубы при положительном расходе
    double rho_out = 850; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе
    
    double visc_in = 15e-6; // плотность нефти, закачиваемой на входе трубы при положительном расходе
    double visc_out = 5e-6; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе



    vector<vector<layer_t>> profiles = vector<vector<layer_t>>(3);

    for (size_t index = 0; index < 100; ++index)
    {

        single_var_moc_t& density_layer_prev = buffer_density.previous();
        single_var_moc_t& density_layer_next = buffer_density.current();

        single_var_moc_t& viscosity_layer_prev = buffer_viscosity.previous();
        single_var_moc_t& viscosity_layer_next = buffer_viscosity.current();

        single_var_moc_t& pressure_layer_prev = buffer_pressure.previous();
        single_var_moc_t& pressure_layer_next = buffer_pressure.current();

        moc_solver<1> solver_density(advection_model, density_layer_prev, density_layer_next);
        moc_solver<1> solver_viscosity(advection_model, viscosity_layer_prev, viscosity_layer_next);

        double dt = solver_density.prepare_step();
        solver_density.step_optional_boundaries(dt, rho_in, rho_out);

        dt = solver_viscosity.prepare_step();
        solver_viscosity.step_optional_boundaries(dt, visc_in, visc_out);

        /// Создаем расчетную модель трубы
        Pipe_model_oil_parties_t pipeModel(pipe, density_layer_next.vars.point_double[0],
            viscosity_layer_next.vars.point_double[0], Q[0]);

        solve_euler_corrector<1>(pipeModel, 1, Pin, &pressure_layer_next.vars.point_double[0]);


        buffer_density.advance(+1);
        buffer_viscosity.advance(+1);
        buffer_pressure.advance(+1);

        // Запись профилей для вывода в файл
        profiles[0].push_back(density_layer_next.vars.point_double[0]);
        profiles[1].push_back(viscosity_layer_next.vars.point_double[0]);
        profiles[2].push_back(pressure_layer_next.vars.point_double[0]);

    }

    // Заголовки для вывода профилей в файл
    vector<string> profile_names = vector<string>();
    profile_names.push_back("плотность [кг/м3]");
    profile_names.push_back("вязкость [Па*с]");
    profile_names.push_back("давление [Па]");

    // Вывод в файл
    profiles_to_file_(pipe, profiles, profile_names);
}
