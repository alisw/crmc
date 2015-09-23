INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/pythia")
FILE(GLOB files *6.4.*.f)
FILE(GLOB files2 *epos.f)
set_source_files_properties(${files} ${files2} ../crmc-aaa.f ../models.F PROPERTIES COMPILE_FLAGS "-D __PYTHIA__")

IF (__CRMCSTATIC__)
add_library(Pythia STATIC ${files} ${files2} ../crmc-aaa.f ../models.F)
ELSE (__CRMCSTATIC__)
add_library(Pythia SHARED ${files} ${files2} ../crmc-aaa.f ../models.F)
target_link_libraries(Pythia CrmcBasic)
ENDIF (__CRMCSTATIC__)


INSTALL (TARGETS Pythia
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
