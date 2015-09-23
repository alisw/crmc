INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/qgsjetII")
FILE(GLOB files *.f ../lzma-read/*.cc)
set_source_files_properties( ${files} ../crmc-aaa.f ../models.F PROPERTIES COMPILE_FLAGS "-D __QGSJETII03__")

IF (__CRMCSTATIC__)
add_library(QgsjetII03 STATIC ${files})
ELSE (__CRMCSTATIC__)
add_library(QgsjetII03 SHARED ${files} ../crmc-aaa.f ../models.F)
target_link_libraries(QgsjetII03 CrmcBasic)
ENDIF (__CRMCSTATIC__)

INSTALL (TARGETS QgsjetII03
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
