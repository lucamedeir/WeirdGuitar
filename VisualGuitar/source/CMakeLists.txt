cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)

project(WeirdGuitar)


add_library(sdl2 STATIC IMPORTED)
set_target_properties(sdl2
	PROPERTIES
	IMPORTED_LOCATION "${LIB_PATH}/SDL2.lib"
	INTERFACE_INCLUDE_DIRECTORIES "${INCLUDE_PATH}"
	)

add_library(gsl STATIC IMPORTED)
set_target_properties(gsl
	PROPERTIES
	IMPORTED_LOCATION "${LIB_PATH}/gsl.lib"
	INTERFACE_INCLUDE_DIRECTORIES "${INCLUDE_PATH}"
	)


add_executable(WeirdGuitar
  main.cpp
)

target_link_libraries(WeirdGuitar sdl2 gsl)

