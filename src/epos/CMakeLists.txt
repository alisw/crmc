INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/epos")
FILE(GLOB files *.f)

IF (__CRMCSTATIC__)
add_library(CrmcBasic STATIC ${files}) #always needed
INSTALL (TARGETS CrmcBasic
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
ELSE (__CRMCSTATIC__)
add_library(CrmcBasic SHARED ${files}) #always needed
INSTALL (TARGETS CrmcBasic
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
add_library(Epos SHARED ../crmc-aaa.f ../models.F)
INSTALL (TARGETS Epos
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
target_link_libraries(Epos CrmcBasic)
ENDIF (__CRMCSTATIC__)

