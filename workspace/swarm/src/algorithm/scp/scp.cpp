// #include <cstdio>
// #include <ament_index_cpp/get_package_share_directory.hpp>
// #include <yaml-cpp/yaml.h>
// #include <string>

// inline std::string get_config_path(const std::string & yaml_name)
// {
//   // this will give you ".../<install_prefix>/share/swarm"
//   //auto pkg_share = rclcpp::get_package_share_directory("swarm");
//   auto pkg_share = ament_index_cpp::get_package_share_directory("swarm");
//   return pkg_share + "/config/" + yaml_name;
// }

// int main(int argc, char ** argv)
// {
//   auto path = get_config_path("example.yaml");
//   YAML::Node cfg = YAML::LoadFile(path);
//   (void) argc;
//   (void) argv;

//   printf("hello world swarm package\n");
//   return 0;
// }
   
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include "yaml-cpp/yaml.h"

#include "scp.hpp"
#include <eigen-quadprog/QuadProg.h>
#include <eigen-quadprog/eigen_quadprog_api.h>



float binomialCoeff(float n, float k){
  if(k==0 || k==n){return 1;}
  return binomialCoeff(n-1,k-1) + binomialCoeff(n-1,k);
}

five_var bernsteinCoeffOrder10(float n, float t_min, float t_max, AXXf t_actual, int num){

  five_var s;
  float l = t_max - t_min;
  AXXf t = (t_actual-t_min)/l;   // scaling the t to (0 to 1)
  AXXf P(num,(int)n+1), Pdot(num,(int)n+1), Pddot(num,(int)n+1); // n is the order of the polynominal , n order has n+1 coefficient
  AXXf Pdddot(num,(int)n+1),Pddddot(num,(int)n+1);
  // P
  for(int i=0;i<=n;++i){
    P.col(i) = binomialCoeff(n,i)*pow(1-t,n-i)*pow(t,i);
  }
  
  // Pdot
  Pdot.col(0) = -10.0 * pow(-t + 1, 9);
  Pdot.col(1) = -90.0 * t * pow(-t + 1, 8) + 10.0 * pow(-t + 1, 9);
  Pdot.col(2) = -360.0 * pow(t, 2) * pow(-t + 1, 7) + 90.0 * t * pow(-t + 1, 8);
  Pdot.col(3) = -840.0 * pow(t, 3) * pow(-t + 1, 6) + 360.0 * pow(t, 2) * pow(-t + 1, 7);
  Pdot.col(4) = -1260.0 * pow(t, 4) * pow(-t + 1, 5) + 840.0 * pow(t, 3) * pow(-t + 1, 6);
  Pdot.col(5) = -1260.0 * pow(t, 5) * pow(-t + 1, 4) + 1260.0 * pow(t, 4) * pow(-t + 1, 5);
  Pdot.col(6) = -840.0 * pow(t, 6) * pow(-t + 1, 3) + 1260.0 * pow(t, 5) * pow(-t + 1, 4);
  Pdot.col(7) = -360.0 * pow(t, 7) * pow(-t + 1, 2) + 840.0 * pow(t, 6) * pow(-t + 1, 3);
  Pdot.col(8) = 45.0 * pow(t, 8) * (2 * t - 2) + 360.0 * pow(t, 7) * pow(-t + 1, 2);
  Pdot.col(9) = -10.0 * pow(t, 9) + 9 * pow(t, 8) * (-10.0 * t + 10.0);
  Pdot.col(10) = 10.0 * pow(t, 9);

  Pddot.col(0) = 90.0 * pow(-t + 1, 8.0);
  Pddot.col(1) = 720.0 * t * pow(-t + 1, 7) - 180.0 * pow(-t + 1, 8);
  Pddot.col(2) = 2520.0 * pow(t, 2) * pow(-t + 1, 6) - 1440.0 * t * pow(-t + 1, 7) + 90.0 * pow(-t + 1, 8);
  Pddot.col(3) = 5040.0 * pow(t, 3) * pow(-t + 1, 5) - 5040.0 * pow(t, 2) * pow(-t + 1, 6) + 720.0 * t * pow(-t + 1, 7);
  Pddot.col(4) = 6300.0 * pow(t, 4) * pow(-t + 1, 4) - 10080.0 * pow(t, 3) * pow(-t + 1, 5) + 2520.0 * pow(t, 2) * pow(-t + 1, 6);
  Pddot.col(5) = 5040.0 * pow(t, 5) * pow(-t + 1, 3) - 12600.0 * pow(t, 4) * pow(-t + 1, 4) + 5040.0 * pow(t, 3) * pow(-t + 1, 5);
  Pddot.col(6) = 2520.0 * pow(t, 6) * pow(-t + 1, 2) - 10080.0 * pow(t, 5) * pow(-t + 1, 3) + 6300.0 * pow(t, 4) * pow(-t + 1, 4);
  Pddot.col(7) = -360.0 * pow(t, 7) * (2 * t - 2) - 5040.0 * pow(t, 6) * pow(-t + 1, 2) + 5040.0 * pow(t, 5) * pow(-t + 1, 3);
  Pddot.col(8) = 90.0 * pow(t, 8) + 720.0 * pow(t, 7) * (2 * t - 2) + 2520.0 * pow(t, 6) * pow(-t + 1, 2);
  Pddot.col(9) = -180.0 * pow(t, 8) + 72 * pow(t, 7) * (-10.0 * t + 10.0);
  Pddot.col(10) = 90.0 * pow(t, 8);

  Pdddot.col(0) = -720.0 * pow(-t + 1, 7);
  Pdddot.col(1) = 2160 * pow(1 - t, 7) - 5040 * pow(1 - t, 6) * t;
  Pdddot.col(2) = -15120 * pow(1-t, 5) * pow(t, 2) + 15120 * pow(1 - t, 6) * t - 2160 * pow(1 - t, 7);
  Pdddot.col(3) = -720 * pow(t-1, 4) * (120 * pow(t, 3) - 108 * pow(t, 2) + 24 * t - 1);
  Pdddot.col(4) = 5040 * pow(t-1, 3) * t * (30 * pow(t, 3) - 36*pow(t, 2) + 12 * t - 1);
  Pdddot.col(5) = -15120 * pow(t-1, 2) * pow(t, 2) * (12 * pow(t, 3) - 18*pow(t, 2) + 8 * t - 1);
  Pdddot.col(6) = 5040 * (t - 1) * pow(t, 3) * (30 * pow(t, 3)-54 * pow(t, 2) + 30 * t - 5);
  Pdddot.col(7) = -720 * pow(t, 7) + 10080 * (1-t) * pow(t, 6) - 45360 * pow(1 - t, 2) * pow(t, 5) + 25200 * pow(1-t, 3) * pow(t, 4) - 2520 * pow(t, 6) * (2 * t - 2);
  Pdddot.col(8) = 2160 * pow(t, 7) - 5040 * (1 - t) * pow(t, 6) + 15120 * pow(1-t, 2) * pow(t, 5) + 5040 * pow(t, 6) * (2 * t - 2);
  Pdddot.col(9) = 504 * (10 - 10 * t) * pow(t, 6) - 2160 * pow(t, 7);
  Pdddot.col(10) = 720.0 * pow(t, 7);

  Pddddot.col(0) = 5040.0 * pow(-t +1, 6);
  Pddddot.col(1) = -4320 * pow(t - 1, 5) * (10 * t - 3)-7200 * pow(t - 1, 6);
  Pddddot.col(2) = 10800 * pow(t - 1, 4) * (15 * pow(t, 2) - 9 * t + 1) + 2160 * pow(t-1, 5) * (30 * t - 9);
  Pddddot.col(3) = -20160 * pow(t-1, 3) * (30 * pow(t, 3) - 36 * pow(t, 2) + 12 * t - 1);
  Pddddot.col(4) = 5040 * pow(t-1, 2) * (210 * pow(t, 4) - 336 * pow(t, 3) + 168 * pow(t, 2) - 28 * t + 1);
  Pddddot.col(5) = -30240 * (t - 1) * t * (42 * pow(t, 4) - 84 * pow(t, 3) + 56 * pow(t, 2) - 14 * t + 1);
  Pddddot.col(6) = 1058400 * pow(t, 6) - 2540160 * pow(t, 5) + 2116800 * pow(t, 4) - 705600 * pow(t, 3) + 75600 * pow(t, 2);
  Pddddot.col(7) = -604800 * pow(t, 6) + 1088640 * pow(t, 5) - 604800 * pow(t, 4) + 100800 * pow(t, 3);
  Pddddot.col(8) = 226800 * pow(t, 6) - 272160 * pow(t, 5) + 75600 * pow(t, 4);
  Pddddot.col(9) = 30240 * pow(t, 5) - 50400 * pow(t, 6);
  Pddddot.col(10) = 5040.0 * pow(t, 6);


  s.a = P;
  s.b = Pdot/l;
  s.c = Pddot/(l*l);
	s.d = Pdddot / (l * l * l);
	s.e = Pddddot / (l * l * l * l);
  return s;
}

five_var computeBernstein(AXXf total_time, float t_fin, int num){

  five_var bernsteinMatrix;
  //bernsteinCoeffOrder10(float n, float t_min, float t_max, Eigen::ArrayXXf t_actual, int num)
  bernsteinMatrix = bernsteinCoeffOrder10(10.0,total_time(0),t_fin,total_time,num);
  return bernsteinMatrix;
}

AXXf block_diag(AXXf arr1, AXXf arr2){
  AXXf temp(arr1.rows()+arr2.rows(),arr1.cols()+arr2.cols());
  temp = 0.0;
  temp.topRows(arr1.rows()).leftCols(arr1.cols())= arr1;
  temp.bottomRows(arr2.rows()).rightCols(arr2.cols()) = arr2;
  return temp;
}


AXXf stack(AXXf arr1, AXXf arr2, char ch){
  // vertical stack
  if (ch=='v') {
    AXXf temp(arr1.rows()+arr2.rows(),arr1.cols());
    temp<<arr1, arr2;
    return temp;
  }
  else if( ch=='h'){
    AXXf temp(arr1.rows(),arr1.cols()+arr2.cols());
    temp<<arr1, arr2;
    return temp;
  }
  else{
    std::cout<<"ERROR";
  }
}
int computeXYZ(probData &prob_data, int VERBOSE){
  // Inequality constraints equation 5,6,7,8
  
  // stacking x,y,z inequality constraints (each axis has less than greater than inequality constraints)
  AXXf A_pos_ineq = block_diag(block_diag(stack(prob_data.P,-prob_data.P,'v'),stack(prob_data.P,-prob_data.P,'v')), stack(prob_data.P,-prob_data.P,'v') );
  AXXf A_vel_ineq = block_diag(block_diag(stack(prob_data.Pdot,-prob_data.Pdot,'v'),stack(prob_data.Pdot,-prob_data.Pdot,'v')), stack(prob_data.Pdot,-prob_data.Pdot,'v') );
  AXXf A_acc_ineq =  block_diag(block_diag(stack(prob_data.Pddot,-prob_data.Pddot,'v'),stack(prob_data.Pddot,-prob_data.Pddot,'v')), stack(prob_data.Pddot,-prob_data.Pddot,'v') );

  AXXf b_x_ineq = stack(prob_data.x_max*AXXf::Ones(prob_data.num,1),-prob_data.x_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b
  AXXf b_y_ineq = stack(prob_data.y_max*AXXf::Ones(prob_data.num,1),-prob_data.y_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b
  AXXf b_z_ineq = stack(prob_data.z_max*AXXf::Ones(prob_data.num,1),-prob_data.z_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b
  AXXf b_pos_ineq = stack(stack(b_x_ineq,b_y_ineq,'v'),b_z_ineq,'v');

  AXXf b_vx_ineq = stack(prob_data.vel_max*AXXf::Ones(prob_data.num,1),-prob_data.vel_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b
  AXXf b_vy_ineq = stack(prob_data.vel_max*AXXf::Ones(prob_data.num,1),-prob_data.vel_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b
  AXXf b_vz_ineq = stack(prob_data.vel_max*AXXf::Ones(prob_data.num,1),-prob_data.vel_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b
  AXXf b_vel_ineq = stack(stack(b_vx_ineq,b_vy_ineq,'v'),b_vz_ineq,'v');

  AXXf b_ax_ineq = stack(prob_data.acc_max*AXXf::Ones(prob_data.num,1),-prob_data.acc_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b
  AXXf b_ay_ineq = stack(prob_data.acc_max*AXXf::Ones(prob_data.num,1),-prob_data.acc_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b
  AXXf b_az_ineq = stack(prob_data.acc_max*AXXf::Ones(prob_data.num,1),-prob_data.acc_min*AXXf::Ones(prob_data.num,1),'v');    // Ax>=b ==> -Ax<=b

  // if you want to include the gravity in z direction 
  // if(prob_data.use_thrust_values){

  //   prob_data.acc_max  = (-(2*prob_data.gravity))

  // }
  AXXf b_acc_ineq = stack(stack(b_ax_ineq,b_ay_ineq,'v'),b_az_ineq,'v');
  AXXf A_slack_ineq, b_slack_ineq;

  // Equality constraints for initial conditions
  prob_data.A_eq =   block_diag(\
                            block_diag(\
                            stack(prob_data.P.row(0),stack(prob_data.Pdot.row(0),prob_data.Pddot.row(0),'v'),'v') ,stack(prob_data.P.row(0),stack(prob_data.Pdot.row(0),prob_data.Pddot.row(0),'v'),'v')),\
                            stack(prob_data.P.row(0),stack(prob_data.Pdot.row(0),prob_data.Pddot.row(0),'v'),'v') \
                            );
  prob_data.b_eq = AXXf(9,1); // 3 * ( pos,vel,acc) 
  prob_data.b_eq << prob_data.x_init,prob_data.vx_init,prob_data.ax_init,
                    prob_data.y_init,prob_data.vy_init,prob_data.ay_init,
                    prob_data.z_init,prob_data.vz_init,prob_data.az_init;

  if(prob_data.num_drone!=0){
    if(prob_data.num_static_obs!=0){
      prob_data.num_static_obs += prob_data.num_drone;
      prob_data.x_static_obs =  stack(prob_data.x_static_obs,prob_data.x_drone,'v');
      prob_data.y_static_obs = stack(prob_data.y_static_obs,prob_data.y_drone,'v');
      prob_data.z_static_obs = stack(prob_data.z_static_obs,prob_data.z_drone,'v');

      prob_data.a_static_obs =  stack(prob_data.a_static_obs,prob_data.a_drone,'v');
      prob_data.b_static_obs = stack(prob_data.b_static_obs,prob_data.b_drone,'v');
      prob_data.c_static_obs = stack(prob_data.c_static_obs,prob_data.c_drone,'v');
    }
    else{
      prob_data.num_static_obs += prob_data.num_drone;
      prob_data.x_static_obs =  prob_data.x_drone;
      prob_data.y_static_obs =  prob_data.y_drone;
      prob_data.z_static_obs =  prob_data.z_drone;

      prob_data.a_static_obs =  prob_data.a_drone;
      prob_data.b_static_obs =  prob_data.b_drone;
      prob_data.c_static_obs =  prob_data.c_drone;
    }
  }
	int tries = 0;
	prob_data.qp_fail = 1;
	prob_data.colliding_step = 1000;
	prob_data.weight_lin_slack = prob_data.weight_lin_slack_og;
	prob_data.weight_quad_slack = prob_data.weight_quad_slack_og;
	
	while(tries < 1 && prob_data.qp_fail){		
		
		Eigen :: ArrayXXf temp_cost, temp_lincost;

		if(prob_data.num_static_obs != 0){
			
			// STACK ALL INEQUALITIES
			if(tries == 0)
      {

				if(!prob_data.on_demand)
					continuousCA(prob_data, VERBOSE);


        prob_data.slack = AXXf::Zero(prob_data.num*prob_data.num_static_obs,prob_data.num*prob_data.num_static_obs);
        prob_data.A_eq = block_diag(prob_data.A_eq,AXXf::Zero(prob_data.num*prob_data.num_static_obs,prob_data.num*prob_data.num_static_obs));
        A_pos_ineq = block_diag(A_pos_ineq,AXXf::Zero(prob_data.num*prob_data.num_static_obs,prob_data.num*prob_data.num_static_obs));
        A_vel_ineq = block_diag(A_vel_ineq,AXXf::Zero(prob_data.num*prob_data.num_static_obs,prob_data.num*prob_data.num_static_obs));
        A_acc_ineq = block_diag(A_acc_ineq,AXXf::Zero(prob_data.num*prob_data.num_static_obs,prob_data.num*prob_data.num_static_obs));

        A_slack_ineq = stack(AXXf::Zero(prob_data.num*prob_data.num_static_obs, 3*prob_data.nvar),-Eigen::MatrixXf::Identity(prob_data.num*prob_data.num_static_obs,prob_data.num*prob_data.num_static_obs),'h');

        
        b_pos_ineq = stack(b_pos_ineq, AXXf::Zero(prob_data.num*prob_data.num_static_obs,1),'v');
        b_vel_ineq = stack(b_vel_ineq, AXXf::Zero(prob_data.num*prob_data.num_static_obs,1),'v');
        b_acc_ineq = stack(b_acc_ineq, AXXf::Zero(prob_data.num*prob_data.num_static_obs,1),'v');

        b_slack_ineq = AXXf::Zero(A_slack_ineq.rows(),1);

        prob_data.b_eq = stack(prob_data.b_eq,AXXf::Zero(A_slack_ineq.rows(),1),'v');
      }
      prob_data.A_ineq = stack(stack(stack(A_pos_ineq,A_vel_ineq,'v'),stack(A_acc_ineq,prob_data.A_coll,'v'),'v'),A_slack_ineq,'v');
      prob_data.b_ineq = stack(stack(stack(b_pos_ineq,b_vel_ineq,'v'),stack(b_acc_ineq,prob_data.b_coll,'v'),'v'),b_slack_ineq,'v');


      temp_cost = block_diag(prob_data.weight_goal*prob_data.cost_goal + prob_data.weight_smoothness*prob_data.cost_smoothness,prob_data.weight_quad_slack*prob_data.cost_slack);
      temp_lincost = stack(-prob_data.weight_goal*prob_data.lincost_goal,prob_data.weight_lin_slack*AXXf::Ones(prob_data.num*prob_data.num_static_obs, 1),'v');
    }
    else{
      prob_data.A_ineq = stack(stack(A_pos_ineq,A_vel_ineq,'v'),A_acc_ineq,'v');
      prob_data.b_ineq = stack(stack(b_pos_ineq,b_vel_ineq,'v'),A_acc_ineq,'v');
      temp_cost = prob_data.weight_goal*prob_data.cost_goal + prob_data.weight_smoothness*prob_data.cost_smoothness;
      temp_lincost = -prob_data.weight_goal*prob_data.lincost_goal;

    }

    // solve qp
    Eigen::QuadProgDense solver_xyz(temp_cost.rows(), prob_data.A_eq.rows(), prob_data.A_ineq.rows());

		solver_xyz.solve((temp_cost).cast<double>().matrix(), temp_lincost.cast<double>().matrix(), 
						prob_data.A_eq.cast<double>().matrix(), prob_data.b_eq.cast<double>().matrix(),
						prob_data.A_ineq.cast<double>().matrix(), prob_data.b_ineq.cast<double>().matrix());
		
    Eigen :: ArrayXXf sol = solver_xyz.result().cast<float>();
		Eigen :: ArrayXXf sol_xy = sol.topRows(2*prob_data.nvar);
		Eigen :: ArrayXXf sol_xyz = sol.topRows(3*prob_data.nvar);

		Eigen :: ArrayXXf sol_x = sol_xy.topRows(prob_data.nvar);
		Eigen :: ArrayXXf sol_y = sol_xy.bottomRows(prob_data.nvar);
		Eigen :: ArrayXXf sol_z = sol_xyz.bottomRows(prob_data.nvar);

    if(prob_data.num_static_obs!=0){
      prob_data.slack = sol.bottomRows(prob_data.num*prob_data.num_static_obs);
    }
		prob_data.x = prob_data.P.matrix() * sol_x.matrix();
		prob_data.y = prob_data.P.matrix() * sol_y.matrix();
		prob_data.z = prob_data.P.matrix() * sol_z.matrix();

		prob_data.xdot = prob_data.Pdot.matrix() * sol_x.matrix();
		prob_data.ydot = prob_data.Pdot.matrix() * sol_y.matrix();
		prob_data.zdot = prob_data.Pdot.matrix() * sol_z.matrix();

		prob_data.xddot = prob_data.Pddot.matrix() * sol_x.matrix();
		prob_data.yddot = prob_data.Pddot.matrix() * sol_y.matrix();
		prob_data.zddot = prob_data.Pddot.matrix() * sol_z.matrix();

		prob_data.xdddot = prob_data.Pdddot.matrix() * sol_x.matrix();
		prob_data.ydddot = prob_data.Pdddot.matrix() * sol_y.matrix();
		prob_data.zdddot = prob_data.Pdddot.matrix() * sol_z.matrix();

		prob_data.xddddot = prob_data.Pddddot.matrix() * sol_x.matrix();
		prob_data.yddddot = prob_data.Pddddot.matrix() * sol_y.matrix();
		prob_data.zddddot = prob_data.Pddddot.matrix() * sol_z.matrix();
		
    prob_data.x_up = prob_data.P_up.matrix() * sol_x.matrix();
		prob_data.y_up = prob_data.P_up.matrix() * sol_y.matrix();
		prob_data.z_up = prob_data.P_up.matrix() * sol_z.matrix();

		prob_data.xdot_up = prob_data.Pdot_up.matrix() * sol_x.matrix();
		prob_data.ydot_up = prob_data.Pdot_up.matrix() * sol_y.matrix();
		prob_data.zdot_up = prob_data.Pdot_up.matrix() * sol_z.matrix();

		prob_data.xddot_up = prob_data.Pddot_up.matrix() * sol_x.matrix();
		prob_data.yddot_up = prob_data.Pddot_up.matrix() * sol_y.matrix();
		prob_data.zddot_up = prob_data.Pddot_up.matrix() * sol_z.matrix();

		prob_data.qp_fail = solver_xyz.fail();
		tries++;
  }

}

AXXf reshape(AXXf x, uint32_t num_r, uint32_t num_c){
  Eigen::Map<Eigen::ArrayXXf> rx(x.data(),num_r,num_c);
  return rx;
}
void continousCA(probData &prob_data, int VERBOSE){

  prob_data.cost_slack  = AXXf(prob_data.num*prob_data.num_static_obs,prob_data.num*prob_data.num_static_obs); // quadratic matrix 0.5(e.TQe) + q.Te
  prob_data.cost_slack.matrix().setIdentity(); // converted to matrix for setting the 2d array to identity
  if(prob_data.world == 2 )
  {

    // 2d space
    AXXf A_f, temp_x, temp_y, A_coll;
    // prob_data.x is a coloumn vector  linearizing the constraints , each row is the linearized obstacle avoidance constraints accross the horizon length.
    // different rows are for the different obstacle
    AXXf temp_delta_fx = -2.0*( (-prob_data.x_static_obs).rowwise() + prob_data.x.transpose().row(0))/pow(prob_data.a_static_obs,2);
    
    // transpose is required, since Eigen is col major by default. we plan to stack the constraints of ego vehicle to the 1st obstacle for the entire horizon. then for the second obstacle and so on  ....
    AXXf temp_A_fx = reshape(temp_delta_fx.transpose(),prob_data.num_static_obs*prob_data.num, 1);
    

    AXXf temp_delta_fy = -2.0*( (-prob_data.y_static_obs).rowwise() + prob_data.y.transpose().row(0))/pow(prob_data.b_static_obs,2);
    AXXf temp_A_fy = reshape(temp_delta_fy.transpose(),prob_data.num_static_obs*prob_data.num, 1);
    for(int i=0;i< prob_data.num_static_obs;++i){
      if(i==0){
        prob_data.A_obs = prob_data.P;
        temp_x = prob_data.x;
        temp_y = prob_data.y;
      }
      else{
        prob_data.A_obs = stack(prob_data.A_obs,prob_data.P,'v');
        temp_x = stack(temp_x, prob_data.x,'v') ;
        temp_y  = stack(temp_y, prob_data.y,'v') ;
      }
    }
    // gradient wrt to the Pcx (prob_data.A_obs.colwise()*temp_A_fx.col(0)) , then stack the y-axis gradient along the colomn
    A_f = stack(prob_data.A_obs.colwise()*temp_A_fx.col(0), prob_data.A_obs.colwise()*temp_A_fy.col(0),'h');  

    A_coll = stack(A_f,-Eigen::MatrixXf::Identity(A_f.rows(),A_f.cols()),'h');

    AXXf temp_b_coll = -(1- ( \
                  pow((-prob_data.x_static_obs).rowwise() + prob_data.x.transpose().row(0),2)/pow(prob_data.a_static_obs, 2) +\
                  pow((-prob_data.y_static_obs).rowwise() + prob_data.y.transpose().row(0),2)/pow(prob_data.b_static_obs, 2) 
                ));
    AXXf  b_coll = reshape(temp_b_coll.transpose(),prob_data.num*prob_data.num_static_obs,1) + temp_A_fx*temp_x + temp_A_fy*temp_y;
    prob_data.A_coll = A_coll;
		prob_data.b_coll = b_coll;
  }
  else{
    // 3d space
    // 2d space
    AXXf A_f, temp_x, temp_y, temp_z, A_coll;
    // prob_data.x is a coloumn vector  linearizing the constraints , each row is the linearized obstacle avoidance constraints accross the horizon length.
    // different rows are for the different obstacle
    AXXf temp_delta_fx = -2.0*( (-prob_data.x_static_obs).rowwise() + prob_data.x.transpose().row(0))/pow(prob_data.a_static_obs,2);
    // transpose is required, since Eigen is col major by default. we plan to stack the constraints of ego vehicle to the 1st obstacle for the entire horizon. then for the second obstacle and so on  ....
    AXXf temp_A_fx = reshape(temp_delta_fx.transpose(),prob_data.num_static_obs*prob_data.num, 1);
    

    AXXf temp_delta_fy = -2.0*( (-prob_data.y_static_obs).rowwise() + prob_data.y.transpose().row(0))/pow(prob_data.b_static_obs,2);
    AXXf temp_A_fy = reshape(temp_delta_fy.transpose(),prob_data.num_static_obs*prob_data.num, 1);

    AXXf temp_delta_fz = -2.0*( (-prob_data.z_static_obs).rowwise() + prob_data.z.transpose().row(0))/pow(prob_data.c_static_obs,2);
    AXXf temp_A_fz = reshape(temp_delta_fz.transpose(),prob_data.num_static_obs*prob_data.num, 1);
    
    for(int i=0;i< prob_data.num_static_obs;++i){
      if(i==0){
        prob_data.A_obs = prob_data.P;
        temp_x = prob_data.x;
        temp_y = prob_data.y;
        temp_z = prob_data.z;

      }
      else{
        prob_data.A_obs = stack(prob_data.A_obs,prob_data.P,'v');
        temp_x = stack(temp_x, prob_data.x,'v') ;
        temp_y  = stack(temp_y, prob_data.y,'v') ;
        temp_z = stack(temp_z, prob_data.z,'v');
      }
    }
    // gradient wrt to the Pcx (prob_data.A_obs.colwise()*temp_A_fx.col(0)) , then stack the y-axis gradient along the colomn
    A_f = stack(stack(prob_data.A_obs.colwise()*temp_A_fx.col(0), prob_data.A_obs.colwise()*temp_A_fy.col(0),'h'),prob_data.A_obs.colwise()*temp_A_fz.col(0),'h');  

    A_coll = stack(A_f,-Eigen::MatrixXf::Identity(A_f.rows(),A_f.cols()),'h');

    AXXf temp_b_coll = -(1- ( \
                  pow((-prob_data.x_static_obs).rowwise() + prob_data.x.transpose().row(0),2)/pow(prob_data.a_static_obs, 2) +\
                  pow((-prob_data.y_static_obs).rowwise() + prob_data.y.transpose().row(0),2)/pow(prob_data.b_static_obs, 2) +\
                  pow((-prob_data.z_static_obs).rowwise() + prob_data.z.transpose().row(0),2)/pow(prob_data.b_static_obs, 2) \
                )
            );
    AXXf  b_coll = reshape(temp_b_coll.transpose(),prob_data.num*prob_data.num_static_obs,1) + temp_A_fx*temp_x + temp_A_fy*temp_y + temp_A_fz*temp_z;
    prob_data.A_coll = A_coll;
		prob_data.b_coll = b_coll;
  }
}