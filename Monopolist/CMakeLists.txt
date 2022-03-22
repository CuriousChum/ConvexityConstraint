cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

project(Monopolist)

set(HeadersDir ${CMAKE_CURRENT_SOURCE_DIR}/Headers)

set(Project_Headers
	${HeadersDir}/Constraint.h
	${HeadersDir}/Constraint.hxx

	${HeadersDir}/ElementaryConstraints.h
	${HeadersDir}/ElementaryConstraints.hxx

	${HeadersDir}/NewtonSolvers.h
	${HeadersDir}/NewtonSolvers.hxx
)

set(Project_Headers_2
	${HeadersDir}/Geometry_2.h
	${HeadersDir}/Geometry_2.hxx
	${HeadersDir}/ElementaryShapes_2.h

	${HeadersDir}/ConvexityConstraint_1.h

	${HeadersDir}/ConvexityConstraint_2.h
	${HeadersDir}/ConvexityConstraint_2.hxx

	${HeadersDir}/LipschitzConstraint_2.h

	${HeadersDir}/PrincipalAgent_2.h
	${HeadersDir}/PrincipalAgent_2.hxx
)

set(Project_Headers_3
	${HeadersDir}/Geometry_3.h
	${HeadersDir}/Geometry_3.hxx
	${HeadersDir}/ElementaryShapes_3.h


	${HeadersDir}/ConvexityConstraint_3.h
	${HeadersDir}/ConvexityConstraint_3.hxx

	${HeadersDir}/PrincipalAgent_3.h
	${HeadersDir}/PrincipalAgent_3.hxx

	${HeadersDir}/Geometry_3.h
	${HeadersDir}/Geometry_3.hxx
)


source_group("Headers" FILES ${Project_Headers})

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${ExternalHeadersDir}) 

add_executable(Monopolist_2 "Monopolist_2_Test.cpp" ${Project_Headers} ${Project_Headers_2})
add_executable(Monopolist_2_NLOpt "Monopolist_2_NLOpt.cpp")
add_executable(Monopolist_3 "Monopolist_3_Test.cpp" ${Project_Headers_3})

target_link_libraries(Monopolist_2 ${GmpLib} ${CgalLib})
target_link_libraries(Monopolist_2_NLOpt ${GmpLib} ${NloptLib} ${CgalLib})
target_link_libraries(Monopolist_3 ${GmpLib} ${NloptLib} ${CgalLib})
