INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/qgsjet")
FILE(GLOB files *.f)
set_source_files_properties( ${files} ../crmc-aaa.f ../models.F PROPERTIES COMPILE_FLAGS "-D __QGSJET01__")

IF (__CRMCSTATIC__)
add_library(Qgsjet01 STATIC ${files})
ELSE (__CRMCSTATIC__)
add_library(Qgsjet01 SHARED ${files} ../crmc-aaa.f ../models.F)
target_link_libraries(Qgsjet01 CrmcBasic)
ENDIF (__CRMCSTATIC__)

INSTALL (TARGETS Qgsjet01
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
