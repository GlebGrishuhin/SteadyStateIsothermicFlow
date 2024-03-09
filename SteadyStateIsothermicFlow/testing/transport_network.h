#pragma once

#include <iomanip>
#include <iostream>
#include <vector>
#include <fstream>
#include <time.h>
#include <algorithm>

#include <fixed/fixed.h>
#include <pde_solvers/pde_solvers.h>

#include <graphlib/graph.h>
#include <graphlib/endogenous.h>
#include <gtest/gtest.h>

using graphlib::edge_t;
using graphlib::graph_t;
using graphlib::net_quickest_ultimate_solver;


namespace transoprt_network
{

    struct task_network_propagation_data_t : graphlib::task_network_propagation_data_t
    {
        static task_network_propagation_data_t default_data()
        {
            const size_t PIPES_NUMBER = 8;
            const size_t T_NUMBER = 3; // номер ребра - отвода на НПЗ

            const vector<double> lengths = { 70.561e3, 79.51e3, 17.93e3, 0.222e3, 59.90e3, 66.79e3, 77.77e3, 56.57e3 };
            const vector<double> diameters = { 1,1,1,0.696,1,1,1,1 };

            task_network_propagation_data_t result;
            // Трубы аналогично топологии ТУ
            result.pipes = vector <pipe_properties_t>();

            simple_pipe_properties simple_LU;

            for (size_t i = 0; i < PIPES_NUMBER; ++i)
            {
                simple_LU = *new simple_pipe_properties();
                simple_LU.length = lengths[i];
                simple_LU.diameter = diameters[i];
                simple_LU.dx = 1000;
                if (i == T_NUMBER)
                    simple_LU.dx = 10;
                result.pipes.push_back(pipe_properties_t::build_simple_pipe(simple_LU));
            }

            vector<double> vol_flows{ 2, 2, 2, 0.1, 1.9, 1.9, 1.9, 1.9 };
            for (size_t index = 0; index < result.pipes.size(); ++index) {
                result.Q.emplace_back(
                    vector<double>(result.pipes[index].profile.getPointCount(), vol_flows[index]));
            }

            for (size_t index = 0; index < result.pipes.size(); ++index) {
                result.models.emplace_back(result.pipes[index], result.Q[index]);
            }

            for (size_t index = 0; index < result.pipes.size(); ++index) {
                result.buffers.emplace_back(2, result.pipes[index].profile.getPointCount());
            }

            result.edges = vector<edge_t>{
                edge_t(0, 1), edge_t(1, 2), edge_t(2, 3),
                edge_t(3, 4), edge_t(3, 5), edge_t(5, 6),
                edge_t(6, 7), edge_t(7, 8) };
            return result;

        }
    };

    class DensityPropagation : public ::testing::Test {
    protected:
        /// @brief Сеть с движением партий для реального трубопровода с разветвлением
        task_network_propagation_data_t net_data;
    protected:
        /// @brief Подготовка к расчету для семейства тестов
        virtual void SetUp() override {
            net_data = task_network_propagation_data_t::default_data();
            double rho_initial = 850;
            net_data.init_buffers(rho_initial);
        }


    };

    /// @brief Расчетный кейс для движения партий с разветвлением
    TEST_F(DensityPropagation, MixDensity) {

        string path = prepare_test_folder();

        std::map<size_t, double> boundaries{
            { 0, 870},
            { 4, -1 },
            { 8, -1 }
        };

        const size_t PRINTED_PIPE = 3; // индекс ребра, для которого будут выводиться профили

        double dt = 300; // шаг по времени

        double T = 300000; // период моделирования
        //double T = 800000; // период моделирования (тест трубы 700км)
        size_t N = static_cast<int>(T / dt);
        double t = 0; // текущее время

        // Временные ряды граничных условий по плотности и расходу
        vector<double> density_in = vector<double>(N, 870);
        vector<double> density_out = vector<double>(N, 870);
        vector<double> volflow_in = vector<double>(N, 2); 
        vector<double> volflow_out = vector<double>(N, 2);



        std::stringstream filename;
        filename << path << "Rho" << ".csv";
        std::ofstream output(filename.str());

        for (size_t index = 0; index < N; ++index) {
            // Учет краевых условий и расходов на новом шаге
            for (size_t i = 1; i < net_data.pipes.size(); ++i) { // Костыль
                double q_pipe = i == 1
                    ? volflow_in[index]
                    : volflow_out[index];
                net_data.Q[i] = vector<double>(net_data.pipes[i].profile.getPointCount(), q_pipe);
            }

            auto& next = net_data.buffers[PRINTED_PIPE].current();
            next.vars.print(t, output);

            boundaries[0] = density_in[index];

            net_quickest_ultimate_solver solver(net_data);
            std::map<size_t, vector<double>> vertices_density
                = solver.step(dt, boundaries);
            net_data.advance_buffers();
            t += dt;
        }
        output.flush();
        output.close();

    }
}

