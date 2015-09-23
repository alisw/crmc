INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/dpmjet/3.0-6")
FILE(GLOB files *.f)
set_source_files_properties( dpmjet3.0-6.f dpmjet-crmc.f phojet1.12-35c4.f pythia6115dpm3v1.f ../../crmc-aaa.f ../../models.F PROPERTIES COMPILE_FLAGS "-D __DPMJET__")

IF (__CRMCSTATIC__)
add_library(Dpmjet STATIC dpmjet3.0-6.f dpmjet-crmc.f phojet1.12-35c4.f pythia6115dpm3v1.f ../../crmc-aaa.f ../../models.F)
ELSE (__CRMCSTATIC__)
add_library(Dpmjet SHARED dpmjet3.0-6.f dpmjet-crmc.f phojet1.12-35c4.f pythia6115dpm3v1.f ../../crmc-aaa.f ../../models.F)
target_link_libraries(Dpmjet CrmcBasic)
ENDIF (__CRMCSTATIC__)

INSTALL (TARGETS Dpmjet
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
