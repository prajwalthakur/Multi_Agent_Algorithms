#!/bin/bash
docker build --rm  $@ -t mp_ros2:latest -f "$(dirname "$0")/../../docker/multi_agent_ros2.Dockerfile" "$(dirname "$0")/../.."