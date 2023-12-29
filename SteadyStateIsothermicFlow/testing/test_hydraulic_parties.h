#pragma once

/// @brief Решение уравнения Бернулли с учетом партий, раcсчитанных методом характеристик
TEST(Hydraulic_Parties, MOC_Euler)
{
    double length = 100;
    double external_diameter = 720e-3;
    double wall_thickness = 10e-3;
    double absolute_roughness = 0.015e-3;
    double Zin = 100;
    double Zout = 50;

    double Pin = 6e6;

    size_t grid_size = 50;

    /// Создаем сущность трубы
    pipe_properties_t pipe;

    /// Задаем сетку трубы
    fill_profile_coordinates_(pipe, length, grid_size);
    fill_profile_heights_(pipe, Zin, Zout);

    pipe.wall.wallThickness = wall_thickness;
    pipe.wall.diameter = external_diameter - 2 * wall_thickness;
    pipe.wall.equivalent_roughness = absolute_roughness;

    typedef profile_collection_t<3> layer_variables_type;
    typedef moc_solver<3>::specific_layer layer_moc_type;

    typedef composite_layer_t<layer_variables_type, layer_moc_type> composite_layer_type;


    typedef profile_collection_t<1> one_param_layer_variables_type;
    typedef moc_solver<1>::specific_layer one_param_layer_moc_type;

    typedef composite_layer_t<one_param_layer_variables_type, one_param_layer_moc_type> one_param_composite_layer_type;

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

    double rho_in = 860; // плотность нефти, закачиваемой на входе трубы при положительном расходе
    double rho_out = 860; // плотность нефти, закачиваемой с выхода трубы при отрицательном расходе 

    double visc_in = 15e-6; // вязкость нефти, закачиваемой на входе трубы при положительном расходе
    double visc_out = 5e-6; // вязкость нефти, закачиваемой с выхода трубы при отрицательном расходе

    // Двумерный вектор граничных условий вида:
    // [[левое ГУ для параметра #1, правое ГУ для параметра #1],
    //  [левое ГУ для параметра #2, правое ГУ для параметра #2]]
    // Номера параметров соответствуют номерам в oil_parameters
    vector<vector<double>> boundary_conditions(NUM_OIL_PROPERTIES, vector<double>(1));
    boundary_conditions[DENSITY_INDEX] = { rho_in, rho_out };
    boundary_conditions[VISCOSITY_INDEX] = { visc_in, visc_out };

    // Ссылки на текущий и предыдущий слои в буфера
    auto& density_layer_prev = std::get<DENSITY_INDEX>(buffer.previous().vars.point_double);
    auto& density_layer_next = std::get<DENSITY_INDEX>(buffer.current().vars.point_double);

    auto& viscosity_layer_prev = std::get<VISCOSITY_INDEX>(buffer.previous().vars.point_double);
    auto& viscosity_layer_next = std::get<VISCOSITY_INDEX>(buffer.current().vars.point_double);

    auto& rho_initial = density_layer_prev;
    rho_initial = vector<double>(rho_initial.size(), 850); // инициализация начальной плотности

    auto& viscosity_initial = viscosity_layer_prev;
    viscosity_initial = vector<double>(viscosity_initial.size(), 10e-6); // инициализация начальной вязкости

    vector<double> Q(pipe.profile.getPointCount(), 10.0); // задаем по трубе расход
    
    // Сущность математической модели адвекции одного параметра
    PipeQAdvection advection_model(pipe, Q);

    // Вектор для хранения исторических профилей и вывода в файл
    vector<vector<layer_t>> profiles = vector<vector<layer_t>>(3);

    for (size_t index = 0; index < 100; ++index)
    {
        // Костыль. Подумать, как убрать. Варианты:
        // 1) Передать в солвер 2 параметра из 3, хранящихся в слое. Как вытащить определенный параметр в виде composite_layer_type?
        // 2) Сделать солвер размерности 3, не обсчитывающий давление. Как? Обнулить коэффициенты нельзя - получим матрицу с det = 0
        // Пробегаем все свойства партии. Для каждого свойства создаем временный слой буфера размерности 1
        // Передаем временный слой в солвер. Пока не ясно
        for (size_t param_index = 0; param_index < oil_parameters.size(); param_index++)
        {
            // Создаем временный слой, куда запишем прошлый профиль данного параметра из буфера
            one_param_composite_layer_type param_layer_prev(pipe.profile.getPointCount());
            param_layer_prev.vars.point_double[0] = buffer.previous().vars.point_double[oil_parameters[param_index]];

            // Создаем временный слой, куда запишем будущий профиль данного параметра из буфера
            one_param_composite_layer_type param_layer_next(pipe.profile.getPointCount());
            param_layer_next.vars.point_double[0] = buffer.current().vars.point_double[oil_parameters[param_index]];

            moc_solver<1> solver(advection_model, param_layer_prev, param_layer_next);
            double dt = solver.prepare_step();

            // Вызов солвера с соответствующими граничными условиями
            solver.step_optional_boundaries(dt, boundary_conditions[oil_parameters[param_index]][0], 
                boundary_conditions[oil_parameters[param_index]][1]);

            // Запись результата расчета в исходный буфер
            buffer.current().vars.point_double[oil_parameters[param_index]] = param_layer_next.vars.point_double[0];
        }
        
       // Мат. модель давления с учетом партийности для передачи в решатель методом Эйлера
       Pipe_model_oil_parties_t pipeModel(pipe, buffer.current().vars.point_double[DENSITY_INDEX],
           buffer.current().vars.point_double[VISCOSITY_INDEX], Q[0]);
       solve_euler_corrector<1>(pipeModel, 1, Pin, &buffer.current().vars.point_double[PRESSURE_INDEX]);
       buffer.advance(+1);

       // Запись профилей для вывода в файл
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
    profiles_to_file_(pipe, profiles, profile_names);

}