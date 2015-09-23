INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/hijing")
FILE(GLOB files hi*.f)
#FILE(GLOB filesPythia ../pythia/*6.4.22.f ../pythia/*epos.f)
FILE(GLOB filesPythia)

set_source_files_properties( ${files} ${filesPythia} ../crmc-aaa.f ../models.F PROPERTIES COMPILE_FLAGS "-D __HIJING__")

IF (__CRMCSTATIC__)
add_library(Hijing STATIC ${files} ${files2} ${filesPythia} ../crmc-aaa.f ../models.F)
ELSE (__CRMCSTATIC__)
add_library(Hijing SHARED ${files} ${files2} ${filesPythia} ../crmc-aaa.f ../models.F)
target_link_libraries(Hijing CrmcBasic)
ENDIF (__CRMCSTATIC__)

INSTALL (TARGETS Hijing
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
