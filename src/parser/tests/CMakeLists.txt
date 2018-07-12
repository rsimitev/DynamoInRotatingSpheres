
# Test parsing of a line.
add_executable(test_parse test_parse.f90)
target_link_libraries(test_parse parser)
add_test(NAME parse COMMAND test_parse
   ${CMAKE_CURRENT_SOURCE_DIR}/data/parse_test.txt)
set_tests_properties(parse PROPERTIES PASS_REGULAR_EXPRESSION "Pass.")

