cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

project(Minkowski)

set(HeadersDir ${CMAKE_CURRENT_SOURCE_DIR}/Headers)
set(Project_Headers
	${HeadersDir}/Meissner_2.h
	${HeadersDir}/Meissner_2.hxx

	${HeadersDir}/Meissner_3.h
	${HeadersDir}/Meissner_3.hxx

	${HeadersDir}/Minkowski_2.h
	${HeadersDir}/Minkowski_2.hxx

	${HeadersDir}/Minkowski_3.h
	${HeadersDir}/Minkowski_3.hxx

	${HeadersDir}/QuotientedNewton.h
	${HeadersDir}/QuotientedNewton.hxx

	${HeadersDir}/ConstraintsProduct.h
	)
source_group("Headers" FILES ${Project_Headers})

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${ExternalHeadersDir}) 

add_executable(Minkowski_2 "Minkowski_2_Test.cpp" ${Project_Headers})
add_executable(Minkowski_3 "Minkowski_3_Test.cpp")
add_executable(Meissner_2 "Meissner_2_Test.cpp")
add_executable(Meissner_3 "Meissner_3_Test.cpp")

target_link_libraries(Meissner_2 ${GmpLib} ${NloptLib})
target_link_libraries(Meissner_3 ${GmpLib} ${NloptLib})
