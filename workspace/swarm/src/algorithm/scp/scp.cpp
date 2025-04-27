#include <cstdio>
#include <ament_index_cpp/get_package_share_directory.hpp>
#include <yaml-cpp/yaml.h>
#include <string>

inline std::string get_config_path(const std::string & yaml_name)
{
  // this will give you ".../<install_prefix>/share/swarm"
  //auto pkg_share = rclcpp::get_package_share_directory("swarm");
  auto pkg_share = ament_index_cpp::get_package_share_directory("swarm");
  return pkg_share + "/config/" + yaml_name;
}

int main(int argc, char ** argv)
{
  auto path = get_config_path("example.yaml");
  YAML::Node cfg = YAML::LoadFile(path);
  (void) argc;
  (void) argv;

  printf("hello world swarm package\n");
  return 0;
}
