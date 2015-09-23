INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/sibyll")
FILE(GLOB files *.f)
set_source_files_properties( ${files} ../crmc-aaa.f ../models.F PROPERTIES COMPILE_FLAGS "-D __SIBYLL__")

IF (__CRMCSTATIC__)
add_library(Sibyll STATIC ${files})
ELSE (__CRMCSTATIC__)
add_library(Sibyll SHARED ${files} ../crmc-aaa.f ../models.F)
target_link_libraries(Sibyll CrmcBasic)
ENDIF (__CRMCSTATIC__)

INSTALL (TARGETS Sibyll
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
