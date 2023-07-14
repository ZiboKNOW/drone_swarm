cmake_minimum_required(VERSION 3.16)
project(Uwb_Location)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  message_generation
  roscpp
  rospy
  serial
  std_msgs
)

find_package(Ceres REQUIRED)

## Generate messages in the 'msg' folder
add_message_files(
  FILES
  uwb.msg
  
)

## Generate added messages and services with any dependencies listed here
generate_messages(
  DEPENDENCIES
  std_msgs
)


## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES Uwb_Location
CATKIN_DEPENDS roscpp rospy std_msgs message_runtime serial
#  DEPENDS system_lib
)



## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
  include
  ${catkin_INCLUDE_DIRS}
)

include_directories(${CERES_INCLUDE_DIRS})

add_library(trilateration include/Uwb_Location/trilateration.h src/trilateration.cpp)
add_dependencies(trilateration ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(trilateration
         ${catkin_LIBRARIES} ${CERES_LIBRARIES})


#可执行文件配置
add_executable(UwbLocation src/UwbLocation.cpp src/trilateration.cpp )
add_dependencies(UwbLocation ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})
target_link_libraries(UwbLocation
        ${catkin_LIBRARIES} ${CERES_LIBRARIES})
