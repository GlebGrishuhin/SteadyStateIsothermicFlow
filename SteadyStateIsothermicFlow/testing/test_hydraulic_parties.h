#pragma once

/// @brief Решение уравнения Бернулли с учетом партий, раcсчитанных методом характеристик
TEST(Hydraulic_Parties, MOC_Euler)
{
    double length = 100000;
    double external_diameter = 720e-3;
    double wall_thickness = 10e-3;
    double absolute_roughness = 0.015e-3;
    double Zin = 100;
    double Zout = 50;

    double Pin = 6e6;

    size_t grid_size = 100;

    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    fill_profile_coordinates_(pipe, length, grid_size);
    fill_profile_heights_(pipe, Zin, Zout);

    pipe.wall.wallThickness = wall_thickness;
    pipe.wall.diameter = external_diameter - 2 * wall_thickness;
    pipe.wall.equivalent_roughness = absolute_roughness;


    typedef composite_layer_t<profile_collection_t<3>, moc_solver<1>::specific_layer> 
        composite_layer_type;

    // Буфер для плотности - индекс [0],
    // вязкости - индекс [1], 
    // давления - индекс [2]
    ring_buffer_t<composite_layer_type> buffer(2, pipe.profile.getPointCount());

    // Константы для навигации по слоям в буфере
    const size_t DENSITY_INDEX = 0;
    const size_t VISCOSITY_INDEX = 1;
    const size_t PRESSURE_INDEX = 2;
    // Количество свойств партии и массив индексов свойств для итерирования по всем свойствам
    const size_t NUM_OIL_PROPERTIES = 2;
    array oil_parameters{ DENSITY_INDEX, VISCOSITY_INDEX };

    double rho_in = 800; // плотность нефти, закачиваемой на входе трубы при положительном расходе
    double rho_out = 860; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе 

    double visc_in = 10e-6; // вязкость нефти, закачиваемой на входе трубы при положительном расходе
    double visc_out = 5e-6; // вязкость нефти, закачиваемой с выхода трубы при отрицательном расходе

    // Двумерный вектор граничных условий вида:
    // [[левое ГУ для параметра #1, правое ГУ для параметра #1],
    //  [левое ГУ для параметра #2, правое ГУ для параметра #2]]
    // Номера параметров соответствуют номерам в oil_parameters
    vector<vector<double>> boundary_conditions(NUM_OIL_PROPERTIES, vector<double>(1));
    boundary_conditions[DENSITY_INDEX] = { rho_in, rho_out };
    boundary_conditions[VISCOSITY_INDEX] = { visc_in, visc_out };

    // Ссылки на текущий и предыдущий слои в буфера
    auto& density_layer_prev = buffer.previous().vars.point_double[DENSITY_INDEX];
    auto& density_layer_next = buffer.current().vars.point_double[DENSITY_INDEX];

    auto& viscosity_layer_prev = std::get<VISCOSITY_INDEX>(buffer.previous().vars.point_double);
    auto& viscosity_layer_next = std::get<VISCOSITY_INDEX>(buffer.current().vars.point_double);

    auto& rho_initial = density_layer_prev;
    rho_initial = vector<double>(rho_initial.size(), 900); // инициализация начальной плотности

    auto& viscosity_initial = viscosity_layer_prev;
    viscosity_initial = vector<double>(viscosity_initial.size(), 15e-6); // инициализация начальной вязкости

    vector<double> Q(pipe.profile.getPointCount(), 0.2); // задаем по трубе расход

    // Сущность математической модели адвекции одного параметра
    PipeQAdvection advection_model(pipe, Q);

    // Вектор для хранения исторических профилей и вывода в файл
    vector<vector<layer_t>> profiles = vector<vector<layer_t>>(3);
    vector<double> time_list;
    time_list.push_back(0);

    for (size_t index = 0; index < 100; ++index)
    {
        auto& prev = buffer.previous();
        auto& next = buffer.current();

        moc_layer_wrapper<1> moc_prev_density(prev.vars.point_double[DENSITY_INDEX],  
            std::get<0>(prev.specific));
        moc_layer_wrapper<1> moc_next_density(next.vars.point_double[DENSITY_INDEX],
            std::get<0>(next.specific));

        moc_solver<1> solver(advection_model, moc_prev_density, moc_next_density);
        double dt = solver.prepare_step();

        //// Вызов солвера с соответствующими граничными условиями
        solver.step_optional_boundaries(dt, boundary_conditions[oil_parameters[DENSITY_INDEX]][0],
            boundary_conditions[oil_parameters[DENSITY_INDEX]][1]);

        moc_layer_wrapper<1> moc_prev_viscosity(prev.vars.point_double[VISCOSITY_INDEX],
            std::get<0>(prev.specific));
        moc_layer_wrapper<1> moc_next_viscosity(next.vars.point_double[VISCOSITY_INDEX],
            std::get<0>(next.specific));

        moc_solver<1> solver2(advection_model, moc_prev_viscosity, moc_next_viscosity);
        dt = solver2.prepare_step();

        //// Вызов солвера с соответствующими граничными условиями
        solver2.step_optional_boundaries(dt, boundary_conditions[oil_parameters[VISCOSITY_INDEX]][0],
            boundary_conditions[oil_parameters[VISCOSITY_INDEX]][1]);


        // Мат. модель давления с учетом партийности для передачи в решатель методом Эйлера
        Pipe_model_oil_parties_t pipeModel(pipe, buffer.current().vars.point_double[DENSITY_INDEX],
            buffer.current().vars.point_double[VISCOSITY_INDEX], Q[0]);
        solve_euler_corrector<1>(pipeModel, 1, Pin, &buffer.current().vars.point_double[PRESSURE_INDEX]);
       
        buffer.advance(+1);

        // Запись профилей для вывода в файл
        time_list.push_back(time_list.back() + dt);
        profiles[0].push_back(buffer.previous().vars.point_double[DENSITY_INDEX]);
        profiles[1].push_back(buffer.previous().vars.point_double[VISCOSITY_INDEX]);
        profiles[2].push_back(buffer.previous().vars.point_double[PRESSURE_INDEX]);
    }

    // Заголовки для вывода профилей в файл
    vector<string> profile_names = vector<string>();
    profile_names.push_back("плотность [кг/м3]");
    profile_names.push_back("вязкость [Па*с]");
    profile_names.push_back("давление [Па]");

    // Вывод в файл
    profiles_to_file_(pipe, time_list, profiles, profile_names);

}