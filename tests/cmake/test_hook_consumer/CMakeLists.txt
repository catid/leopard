cmake_minimum_required(VERSION 3.7)
project(leopard_test_hook_consumer LANGUAGES CXX)

find_package(leopard CONFIG REQUIRED)
if(NOT TARGET leopard::leopard)
    message(FATAL_ERROR "Installed package did not export leopard::leopard")
endif()

add_executable(leopard_test_hook_consumer main.cpp)
target_link_libraries(leopard_test_hook_consumer PRIVATE leopard::leopard)
