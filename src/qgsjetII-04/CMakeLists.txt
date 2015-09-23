INCLUDE_DIRECTORIES ("${PROJECT_SOURCE_DIR}/src/qgsjetII-04")
FILE(GLOB files *.f ../lzma-read/*.cc)
set_source_files_properties( ${files} ../crmc-aaa.f ../models.F PROPERTIES COMPILE_FLAGS "-D __QGSJETII04__")

IF (__CRMCSTATIC__)
add_library(QgsjetII04 STATIC ${files})
ELSE (__CRMCSTATIC__)
add_library(QgsjetII04 SHARED ${files} ../crmc-aaa.f ../models.F)
target_link_libraries(QgsjetII04 CrmcBasic)
ENDIF (__CRMCSTATIC__)


INSTALL (TARGETS QgsjetII04
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
        )
