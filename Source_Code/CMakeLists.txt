cmake_minimum_required(VERSION 3.6)
project(File_Example)

set(CMAKE_C_STANDARD 99)

set(SOURCE_FILES main.c process.c classifiers.c  collection.c training_file.c process_testing_data.c do_not_use_file.c)

add_executable(File_Example ${SOURCE_FILES})
