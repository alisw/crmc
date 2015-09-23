INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/phojet")
FILE(GLOB files *.f)
FILE(GLOB filesPythia ../pythia/pythia6115.f)
set_source_files_properties( ${files} ${filesPythia} ../crmc-aaa.f ../models.F PROPERTIES COMPILE_FLAGS "-D __PHOJET__")

IF (__CRMCSTATIC__)
add_library(Phojet STATIC ${files} ${filesPythia} ../crmc-aaa.f ../models.F)
ELSE (__CRMCSTATIC__)
add_library(Phojet SHARED ${files} ${filesPythia} ../crmc-aaa.f ../models.F)
target_link_libraries(Phojet CrmcBasic)
ENDIF (__CRMCSTATIC__)

INSTALL (TARGETS Phojet
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
