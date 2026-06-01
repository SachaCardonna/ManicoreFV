ADD_EXECUTABLE(debugFV_printmesh printmesh.cpp)
SET_TARGET_PROPERTIES(debugFV_printmesh PROPERTIES OUTPUT_NAME "printMesh")
TARGET_LINK_LIBRARIES(debugFV_printmesh FVmesh_builder quadrature)
ADD_CUSTOM_COMMAND(TARGET debugFV_printmesh POST_BUILD COMMAND 
  ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/printmesh.py ${CMAKE_BINARY_DIR}/printmesh.py)
