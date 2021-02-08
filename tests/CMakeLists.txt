mfront_behaviours_check_library(BehaviourTest
  Elasticity
  Plasticity
  Mazars
  SaintVenantKirchhoffElasticity)

add_executable(PeriodicTestNew
  EXCLUDE_FROM_ALL
  PeriodicTestNew.cxx)
target_link_libraries(PeriodicTestNew
  PRIVATE MFEMMGIS)
  add_dependencies(check PeriodicTestNew)

function(add_periodic_test ncase nsolver)
  set(test "PeriodicTestNew-${ncase}-${nsolver}")
  add_test(NAME ${test}
   COMMAND PeriodicTestNew 
   "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube_2mat_per.mesh"
   "--library" "$<TARGET_FILE:BehaviourTest>"
   "--test-case" "${ncase}"
   "--linearsolver" "${nsolver}")
  if((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest
      PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
  else((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest)
  endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
endfunction(add_periodic_test)

foreach(ncase RANGE 0 5)
  foreach(nsolver RANGE 0 1)
    add_periodic_test(${ncase} ${nsolver})
  endforeach(nsolver)
endforeach(ncase)

add_executable(UnilateralTensileTest
  EXCLUDE_FROM_ALL
  UnilateralTensileTest.cxx)
target_link_libraries(UnilateralTensileTest
  PRIVATE MFEMMGIS)
add_dependencies(check UnilateralTensileTest)

function(add_unilateral_test behaviour internal_state_variable)
  set(test "UnilateralTensileTest-${behaviour}")
  add_test(NAME ${test}
   COMMAND UnilateralTensileTest
   "--mesh" "${CMAKE_CURRENT_SOURCE_DIR}/cube.mesh"
   "--library" "$<TARGET_FILE:BehaviourTest>"
   "--behaviour" "${behaviour}"
   "--reference-file" "${CMAKE_CURRENT_SOURCE_DIR}/references/${behaviour}.ref"
   "--internal-state-variable" "${internal_state_variable}")
  if((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest
      PROPERTY ENVIRONMENT "PATH=$<TARGET_FILE_DIR:MFEMMGIS>\;${MGIS_PATH_STRING}")
  else((CMAKE_HOST_WIN32) AND (NOT MSYS))
    set_property(TEST ${test}
      PROPERTY DEPENDS BehaviourTest)
  endif((CMAKE_HOST_WIN32) AND (NOT MSYS))
endfunction(add_unilateral_test)

add_unilateral_test(Plasticity EquivalentPlasticStrain)
add_unilateral_test(Mazars Damage)
add_unilateral_test(SaintVenantKirchhoffElasticity EquivalentStrain)
