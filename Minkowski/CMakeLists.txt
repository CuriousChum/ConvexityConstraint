cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

project(Minkowski)

set(Project_Headers
	${CMAKE_CURRENT_SOURCE_DIR}/Meissner_2_Test.h
	${CMAKE_CURRENT_SOURCE_DIR}/Meissner_2.h
	${CMAKE_CURRENT_SOURCE_DIR}/Meissner_2.hxx

	${CMAKE_CURRENT_SOURCE_DIR}/Meissner_3_Test.h
	${CMAKE_CURRENT_SOURCE_DIR}/Meissner_3.h
	${CMAKE_CURRENT_SOURCE_DIR}/Meissner_3.hxx

	${CMAKE_CURRENT_SOURCE_DIR}/Minkowski_2_Test.h
	${CMAKE_CURRENT_SOURCE_DIR}/Minkowski_2.h
	${CMAKE_CURRENT_SOURCE_DIR}/Minkowski_2.hxx

	${CMAKE_CURRENT_SOURCE_DIR}/Minkowski_3_Test.h
	${CMAKE_CURRENT_SOURCE_DIR}/Minkowski_3.h
	${CMAKE_CURRENT_SOURCE_DIR}/Minkowski_3.hxx

	${CMAKE_CURRENT_SOURCE_DIR}/QuotientedNewton.h
	${CMAKE_CURRENT_SOURCE_DIR}/QuotientedNewton.hxx
	)
source_group("Headers" FILES ${Project_Headers})

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${ExternalHeadersDir}) 
add_executable(Minkowski "Main.cpp" ${Project_Headers})
target_link_libraries(Minkowski ${GmpLib} ${NloptLib})
