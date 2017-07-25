# TODO add tests
MACRO(ADD_MINC_TEST name command )
 ADD_TEST( ${name}
      ${CMAKE_CURRENT_SOURCE_DIR}/${command} ${CMAKE_CURRENT_BINARY_DIR}/../src ${CMAKE_CURRENT_BINARY_DIR} ${ARGV2} ${ARGV3} ${ARGV4} ${ARGV5} ${ARGV6} ${ARGV7} ${ARGV8} ${ARGV9} ${ARGV10} ${ARGV11}
      ${ARGV12} ${ARGV13} ${ARGV14} ${ARGV15} ${ARGV16} ${ARGV17} ${ARGV18} ${ARGV19} ${ARGV20} ${ARGV21} ${ARGV22}
      ${ARGV23} ${ARGV24} ${ARGV25} ${ARGV26}
 )

 IF(MINC_TEST_ENVIRONMENT)
  set_tests_properties( ${name} PROPERTIES ENVIRONMENT "${MINC_TEST_ENVIRONMENT}")
 ENDIF(MINC_TEST_ENVIRONMENT)
ENDMACRO(ADD_MINC_TEST)

ADD_MINC_TEST(test_nlm  test_nlm.sh )
ADD_MINC_TEST(test_anlm test_anlm.sh )
ADD_MINC_TEST(test_patch_segmentation test_segmentation.sh )
