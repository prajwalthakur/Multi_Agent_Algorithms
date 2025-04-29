#pragma once
#include <eigen3/Eigen/Dense>
#include <vector>
#include "yaml-cpp/yaml.h"
typedef Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> MapMatrixCol;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> MatrixXdRow;
using AXXf = Eigen::ArrayXXf;
using AXf = Eigen::ArrayXf;       // column vector
struct five_var
{
    AXXf a, b, c, d, e;
};       

struct probData
{
    bool mpc, free_space, on_demand;
    int num, num_up, num_ctrl, nvar, kappa, VERBOSE, mpc_step, max_time, num_static_obs, num_drone, qp_fail, id_badge, colliding_step, world;

    float dt, t_plan, weight_smoothness, weight_goal, weight_quad_slack, weight_lin_slack, dist_to_goal, dist_stop, buffer, weight_lin_slack_og, weight_quad_slack_og;
    float vel_max,vel_min, acc_max, acc_min, prox_obs, prox_agent, lx_drone, ly_drone, lz_drone, rmin;

    float x_min, y_min, z_min,
            x_max, y_max, z_max,
            x_init, y_init, z_init,
            x_goal, y_goal, z_goal, 
            vx_init, vy_init, vz_init,
            ax_init, ay_init, az_init;

    float gravity, f_min, f_max;
    float mean, stdev;
    bool use_thrust_values, use_model;

    AXXf A_coll, b_coll;
    AXXf agents_x, agents_y, agents_z; 
    AXXf A_eq, A_ineq, A_obs, cost_goal, cost_smoothness, cost_slack;
    AXXf b_eq, b_ineq, lincost_goal;

    AXXf P, Pdot, Pddot, Pdddot, Pddddot;  
    AXXf P_up, Pdot_up, Pddot_up, Pdddot_up, Pddddot_up;

    AXXf x_ref, y_ref, z_ref;
    AXXf slack;

    AXXf x_static_obs, y_static_obs, z_static_obs;
    AXXf x_static_obs_og, y_static_obs_og, z_static_obs_og;
    AXXf x_drone, y_drone, z_drone;

    AXXf a_drone, b_drone, c_drone;
    AXXf a_static_obs, b_static_obs, c_static_obs;
    AXXf a_static_obs_og, b_static_obs_og, c_static_obs_og; 

    AXXf x, y, z,
        xdot, ydot, zdot,
        xddot, yddot, zddot,
        xdddot, ydddot, zdddot,
        xddddot, yddddot, zddddot;
    
    AXXf x_up, y_up, z_up,
        xdot_up, ydot_up, zdot_up,
        xddot_up, yddot_up, zddot_up;
    
    std :: vector<float> smoothness, arc_length, inter_agent_dist, agent_obs_dist, inter_agent_dist_min, agent_obs_dist_min;
    std :: string solver;
    YAML :: Node params;

    std :: vector<std :: vector<float>> pos_static_obs, dim_static_obs;
};

five_var bernsteinCoeffOrder10(float n, float tmin, float tmax, AXXf t_actual, int num);
five_var computeBernstein(AXXf tot_time, float t_fin, int num);
AXXf stack(AXXf arr1, AXXf arr2, char ch);
AXXf reshape(AXXf x, uint32_t r, uint32_t c);
AXXf arctan2(AXXf arr1, AXXf arr2);
AXXf maximum(float val, AXXf arr2);
AXXf delete_values(float val, AXXf arr);
AXXf block_diag(AXXf arr1, AXXf arr2);
AXXf diff(AXXf arr);

int computeXYZ(probData &prob_data, int VERBOSE);

void initObstacles(probData &prob_data, int VERBOSE);
void neigbhoringAgents(probData &prob_data, int VERBOSE);

void continuousCA(probData &prob_data, int VERBOSE);
void ondemandCA(probData &prob_data, int VERBOSE);
int collidingStep(probData &prob_data, int VERBOSE);

void initializeOptimizer(probData &prob_data, int VERBOSE);
void deployAgent(probData &prob_data, int VERBOSE);
