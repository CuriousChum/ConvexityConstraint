cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)

include(${CMAKE_CURRENT_SOURCE_DIR}/../common.cmake)

project(Monopolist)

#set(HeadersDir ${CMAKE_CURRENT_SOURCE_DIR}/Headers)
set(HeadersDir ${CMAKE_CURRENT_SOURCE_DIR})

set(Project_Headers
	${HeadersDir}/Constraint.h
	${HeadersDir}/Constraint.hxx

	${HeadersDir}/ConvexityConstraint_3.h
	${HeadersDir}/ConvexityConstraint_3.hxx

	${HeadersDir}/ConvexityConstraint.h
	${HeadersDir}/ConvexityConstraint.hxx

#	${HeadersDir}/Minkowski_3.h
#	${HeadersDir}/Minkowski_3.hxx

	${HeadersDir}/Geometry_2.h
	${HeadersDir}/Geometry_2.hxx

	${HeadersDir}/Geometry_3.h
	${HeadersDir}/Geometry_3.hxx

	${HeadersDir}/LipschitzConstraint.h

	${HeadersDir}/NewtonSolvers.h
	${HeadersDir}/NewtonSolvers.hxx

	${HeadersDir}/PrincipalAgent_3.h
	${HeadersDir}/PrincipalAgent_3.hxx

	${HeadersDir}/PrincipalAgent_NLOpt.h
	${HeadersDir}/PrincipalAgent_Test.h

	${HeadersDir}/PrincipalAgent.h
	${HeadersDir}/PrincipalAgent.hxx
	)
source_group("Headers" FILES ${Project_Headers})

include_directories(${PROJECT_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${ExternalHeadersDir}) 

add_executable(Monopolist_2 "main.cpp" ${Project_Headers})
target_link_libraries(Monopolist_2 ${GmpLib} ${NloptLib} ${CgalLib})
