cmake_minimum_required(VERSION 3.12)
project(scloudplus LANGUAGES C)

# 设置项目版本
set(PROJECT_VERSION 1.1.0)

add_compile_options(-O3)
message(STATUS "Enabling O3 optimization")

# 配置选项
option(BUILD_WITH_CLANG "Use Clang compiler" OFF)
option(BUILD_WITH_GCC "Use GCC compiler" ON)

# 默认使用Release构建类型（仅适用于单配置生成器）
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type" FORCE)
    # 设置可用的构建类型选项
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# 设置源文件和头文件
set(SRC_FILES
        src/aes.c
        src/encode.c
        src/random.c
        src/hash.c
        src/sample.c
        src/pke.c
        src/matrix.c
        src/util.c
        src/kem.c
)

set(HEADER_FILES
        include/aes.h
        include/encode.h
        include/random.h
        include/hash.h
        include/sample.h
        include/pke.h
        include/matrix.h
        include/util.h
        include/kem.h
        include/ds_benchmark.h
)

# 检测目标架构 添加 AES 指令集支持
if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    # 添加 AES-NI 指令集支持
    set(AES_FLAGS "-maes")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm")
    # ARM 平台使用 ARM Crypto 扩展
    set(AES_FLAGS "-march=armv8-a+crypto")
else()
    set(AES_FLAGS "")
endif()

# 更新最终标志（保留原有标志，添加AES指令集支持）
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${AES_FLAGS}"
        CACHE STRING "C compiler flags" FORCE)

# 创建三个不同配置的可执行文件
foreach(LEVEL 128 192 256)
    add_executable(scloudplus${LEVEL}_aes
            ${SRC_FILES}
            test/test.c
    )

    target_include_directories(scloudplus${LEVEL}_aes PRIVATE include)
    target_compile_definitions(scloudplus${LEVEL}_aes PRIVATE scloudplus_l=${LEVEL})

endforeach()

# 添加清理目标
add_custom_target(clean-all
        COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}
        COMMENT "Cleaning all build files"
)
