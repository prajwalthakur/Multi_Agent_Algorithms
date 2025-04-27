cmake_minimum_required(VERSION 3.8)
project(swarm)

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(std_msgs REQUIRED)

# find_package(eigen3_cmake_module REQUIRED)
find_package(eigen-quadprog REQUIRED)
find_package(Eigen3            REQUIRED)
find_package(yaml-cpp          REQUIRED)
find_package(eigen-quadprog    REQUIRED)


set(EIGEN_QUADPROG_INCLUDE_DIR "/usr/include/eigen-quadprog")
# include_directories(
#   include

# )

#am
file(GLOB SWARM_AM_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm/am/*.cpp"
)
add_library(lib_am_swarm 
            ${SWARM_AM_SOURCES}
            )  

target_include_directories(lib_am_swarm PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/swarm/am>
            $<INSTALL_INTERFACE:include/swarm/am>
            ${EIGEN_QUADPROG_INCLUDE_DIR}
          )            


ament_target_dependencies(lib_am_swarm
            rclcpp
            yaml-cpp
            
          )           
target_link_libraries(lib_am_swarm  Eigen3::Eigen )



add_executable(swarm_am_nav src/main_am_swarm.cpp)
target_link_libraries(swarm_am_nav lib_am_swarm ) #eigen-quadprog yaml-cpp ${catkin_LIBRARIES})
ament_target_dependencies(swarm_am_nav
  rclcpp
  std_msgs
)

#scp
file(GLOB SWARM_SCP_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/algorithm/scp/*.cpp"
)
add_library(lib_scp_swarm 
            ${SWARM_SCP_SOURCES}
            )   

target_include_directories(lib_scp_swarm PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/swarm/scp>
            $<INSTALL_INTERFACE:include/swarm/scp>
            ${EIGEN_QUADPROG_INCLUDE_DIR}
          )
ament_target_dependencies(lib_scp_swarm
            rclcpp
            yaml-cpp
          )

target_link_libraries(lib_am_swarm  Eigen3::Eigen )          

add_executable(swarm_scp_nav src/main_scp_swarm.cpp)
target_link_libraries(swarm_scp_nav lib_scp_swarm)
ament_target_dependencies(swarm_scp_nav
  rclcpp
  std_msgs
)

install(
  TARGETS
    lib_am_swarm
    swarm_am_nav
    lib_scp_swarm
    swarm_scp_nav
  ARCHIVE DESTINATION lib/${PROJECT_NAME}
  LIBRARY DESTINATION lib/${PROJECT_NAME}
  RUNTIME DESTINATION lib/${PROJECT_NAME}
)


install(
  DIRECTORY include/swarm/am
  DESTINATION include/swarm/
)

install(
  DIRECTORY include/swarm/scp
  DESTINATION include/swarm/
)


install(
  DIRECTORY launch
  DESTINATION share/${PROJECT_NAME}/
)

# install(
#   DIRECTORY 
#     include/swarm/am
#     include/swarm/scp
#   DESTINATION include/swarm
# )

# ament_export_include_directories(
#   include/swarm/am
#   include/swarm/scp
# )


if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # comment the line when a copyright and license is added to all source files
  set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # comment the line when this package is in a git repo and when
  # a copyright and license is added to all source files
  set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

ament_package()
