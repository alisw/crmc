INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/dpmjet/3.2017-1")
FILE(GLOB files *.f)
set_source_files_properties( dpmjet.f dpmjet_epos.f phojet.f pythia-6427_cors.f build.f ../../crmc-aaa.f ../../models.F PROPERTIES COMPILE_FLAGS "-D __DPMJET__")

IF (__CRMCSTATIC__)
add_library(Dpmjet STATIC dpmjet.f dpmjet_epos.f phojet.f pythia-6427_cors.f build.f ../../crmc-aaa.f ../../models.F)
ELSE (__CRMCSTATIC__)
add_library(Dpmjet SHARED dpmjet.f dpmjet_epos.f phojet.f pythia-6427_cors.f build.f ../../crmc-aaa.f ../../models.F)
target_link_libraries(Dpmjet CrmcBasic)
ENDIF (__CRMCSTATIC__)

INSTALL (TARGETS Dpmjet
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
