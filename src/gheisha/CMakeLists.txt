INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/gheisha")
FILE(GLOB files *.f)
set_source_files_properties( ${files} ../crmc-aaa.f ../models.F PROPERTIES COMPILE_FLAGS "-D __GHEISHA__")

IF (__CRMCSTATIC__)
add_library(Gheisha STATIC ${files})
ELSE (__CRMCSTATIC__)
add_library(Gheisha SHARED ${files} ../crmc-aaa.f ../models.F)
target_link_libraries(Gheisha CrmcBasic)
ENDIF (__CRMCSTATIC__)

INSTALL (TARGETS Gheisha
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
