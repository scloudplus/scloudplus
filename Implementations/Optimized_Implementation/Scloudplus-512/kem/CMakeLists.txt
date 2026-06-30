cmake_minimum_required(VERSION 3.11.4)
project(scloudplus_optimized_kem C)

# This CMakeLists.txt is intentionally self-contained for the submitted KEM
# entry directory.  It does not include an external .cmake file: source
# selection, compile definitions, and public test/KAT/benchmark targets are all
# visible from this file.  The actual algorithm parameters are fixed in the
# checked-in parameters.h next to this CMakeLists.txt.

set(SCLOUDPLUS_TIER "Optimized")
get_filename_component(SCLOUDPLUS_LEVEL_PATH "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
get_filename_component(SCLOUDPLUS_LEVEL_DIR "${SCLOUDPLUS_LEVEL_PATH}" NAME)
string(REGEX REPLACE "^Scloudplus-" "" SCLOUDPLUS_LEVEL "${SCLOUDPLUS_LEVEL_DIR}")
set(SCLOUDPLUS_IMPLEMENTATIONS_DIR "${CMAKE_CURRENT_LIST_DIR}/../../..")
set(SCLOUDPLUS_SHARED_DIR "${SCLOUDPLUS_IMPLEMENTATIONS_DIR}/_shared")
set(SCLOUDPLUS_CORE_DIR "${SCLOUDPLUS_SHARED_DIR}/scloudplus_core")
set(SCLOUDPLUS_API_PKC_DIR "${SCLOUDPLUS_SHARED_DIR}/api_pkc")
set(SCLOUDPLUS_TEST_DIR "${SCLOUDPLUS_IMPLEMENTATIONS_DIR}/../Self_Evaluation/tests")
set(SCLOUDPLUS_BENCH_DIR "${SCLOUDPLUS_IMPLEMENTATIONS_DIR}/../Self_Evaluation/benchmark")

option(SCLOUDPLUS_ENABLE_TSV110_TUNE "Enable -mtune=tsv110 for validated Kunpeng AArch64 GNU builds" ON)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(SCLOUDPLUS_FAMILY "AES" CACHE STRING "Scloud+ primitive family: AES, SHAKE, or SM3")
set_property(CACHE SCLOUDPLUS_FAMILY PROPERTY STRINGS AES SHAKE SM3)

set(SCLOUDPLUS_BACKEND "AUTO" CACHE STRING "Scloud+ backend: AUTO, AVX2, or NEON")
set_property(CACHE SCLOUDPLUS_BACKEND PROPERTY STRINGS AUTO AVX2 NEON)

string(TOUPPER "${SCLOUDPLUS_FAMILY}" SCLOUDPLUS_FAMILY)
string(TOUPPER "${SCLOUDPLUS_BACKEND}" SCLOUDPLUS_BACKEND)

set(_valid_family AES SHAKE SM3)
if(NOT SCLOUDPLUS_FAMILY IN_LIST _valid_family)
  message(FATAL_ERROR "SCLOUDPLUS_FAMILY must be AES, SHAKE, or SM3")
endif()

set(_valid_backend AUTO AVX2 NEON)
if(NOT SCLOUDPLUS_BACKEND IN_LIST _valid_backend)
  message(FATAL_ERROR "Optimized builds support SCLOUDPLUS_BACKEND=AUTO, AVX2, or NEON")
endif()
if(SCLOUDPLUS_BACKEND STREQUAL "AUTO")
  if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    set(SCLOUDPLUS_BACKEND AVX2)
  elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
    set(SCLOUDPLUS_BACKEND NEON)
  else()
    message(FATAL_ERROR "Cannot auto-select optimized backend for CMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR}; set SCLOUDPLUS_BACKEND explicitly")
  endif()
  message(STATUS "Auto-selected SCLOUDPLUS_BACKEND=${SCLOUDPLUS_BACKEND} for ${CMAKE_SYSTEM_PROCESSOR}")
endif()


set(SCLOUDPLUS_CORE_SOURCES "")
set(SCLOUDPLUS_KAT_CORE_SOURCES "")
set(SCLOUDPLUS_EXTRA_COMPILE_OPTIONS "")
set(SCLOUDPLUS_EXTRA_COMPILE_DEFINITIONS
  SCLOUDPLUS_FAMILY_${SCLOUDPLUS_FAMILY}
)

function(scloudplus_add_core_source relpath)
  list(APPEND SCLOUDPLUS_CORE_SOURCES "${SCLOUDPLUS_CORE_DIR}/${relpath}")
  if(NOT relpath STREQUAL "common/random.c")
    list(APPEND SCLOUDPLUS_KAT_CORE_SOURCES "${SCLOUDPLUS_CORE_DIR}/${relpath}")
  endif()
  set(SCLOUDPLUS_CORE_SOURCES "${SCLOUDPLUS_CORE_SOURCES}" PARENT_SCOPE)
  set(SCLOUDPLUS_KAT_CORE_SOURCES "${SCLOUDPLUS_KAT_CORE_SOURCES}" PARENT_SCOPE)
endfunction()

set(SCLOUDPLUS_PACK_SOURCE ref/pack_reference.c)

if(SCLOUDPLUS_BACKEND STREQUAL "AVX2")
  if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64")
    message(FATAL_ERROR "AVX2 backend requires an x86_64/AMD64 processor")
  endif()
  list(APPEND SCLOUDPLUS_EXTRA_COMPILE_OPTIONS -O3 -mavx2 -mbmi2 -mpopcnt -maes)
  list(APPEND SCLOUDPLUS_EXTRA_COMPILE_DEFINITIONS
    SCLOUDPLUS_TIER_OPTIMIZED
    SCLOUDPLUS_BACKEND_AVX2
    SCLOUDPLUS_AVX2_FAMILY_${SCLOUDPLUS_FAMILY}
  )
  if(SCLOUDPLUS_FAMILY STREQUAL "SM3")
    list(APPEND SCLOUDPLUS_EXTRA_COMPILE_OPTIONS -funroll-loops -flto)
    set(SCLOUDPLUS_HASH_SOURCE avx2/hash_sm3_avx2.c)
    scloudplus_add_core_source(avx2/sm3_avx2.c)
  else()
    set(SCLOUDPLUS_HASH_SOURCE common/hash_aes_shake.c)
    if(SCLOUDPLUS_FAMILY STREQUAL "AES")
      scloudplus_add_core_source(avx2/aes_aesni.c)
    elseif(SCLOUDPLUS_FAMILY STREQUAL "SHAKE")
      scloudplus_add_core_source(avx2/fips202x4.c)
      scloudplus_add_core_source(avx2/keccak4x/KeccakP-1600-times4-SIMD256.c)
    endif()
  endif()
  set(SCLOUDPLUS_MATRIX_SOURCE avx2/matrix_avx2.c)
  scloudplus_add_core_source(avx2/sample_avx2.c)
elseif(SCLOUDPLUS_BACKEND STREQUAL "NEON")
  if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
    message(FATAL_ERROR "NEON backend requires an AArch64/ARM64 processor")
  endif()
  list(APPEND SCLOUDPLUS_EXTRA_COMPILE_DEFINITIONS
    SCLOUDPLUS_TIER_OPTIMIZED
    SCLOUDPLUS_BACKEND_NEON
    SCLOUDPLUS_NEON_FAMILY_${SCLOUDPLUS_FAMILY}
  )
  if(SCLOUDPLUS_FAMILY STREQUAL "AES")
    list(APPEND SCLOUDPLUS_EXTRA_COMPILE_OPTIONS -O3 -march=armv8-a+crypto)
    set(SCLOUDPLUS_HASH_SOURCE common/hash_aes_shake.c)
    scloudplus_add_core_source(neon/aes_armcrypto.c)
  elseif(SCLOUDPLUS_FAMILY STREQUAL "SHAKE")
    list(APPEND SCLOUDPLUS_EXTRA_COMPILE_OPTIONS -O3 -march=armv8-a+simd)
    set(SCLOUDPLUS_HASH_SOURCE common/hash_aes_shake.c)
  else()
    set(SCLOUDPLUS_HASH_SOURCE common/hash_sm3_portable.c)
    if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
      list(APPEND SCLOUDPLUS_EXTRA_COMPILE_OPTIONS -O3 -march=armv8-a+simd)
    else()
      list(APPEND SCLOUDPLUS_EXTRA_COMPILE_OPTIONS -O3 -march=armv8.2-a+sm4)
    endif()
    scloudplus_add_core_source(neon/sm3_aarch64.c)
  endif()
  set(SCLOUDPLUS_MATRIX_SOURCE neon/matrix_neon.c)
  set(SCLOUDPLUS_PACK_SOURCE neon/pack_neon.c)
  scloudplus_add_core_source(neon/sample_neon.c)
endif()

foreach(relpath
    common/encode.c
    ${SCLOUDPLUS_HASH_SOURCE}
    common/kem.c
    ${SCLOUDPLUS_MATRIX_SOURCE}
    ${SCLOUDPLUS_PACK_SOURCE}
    common/pke.c
    common/random.c
    common/sample.c
    common/util.c)
  scloudplus_add_core_source(${relpath})
endforeach()

if(SCLOUDPLUS_ENABLE_TSV110_TUNE AND
   CMAKE_C_COMPILER_ID STREQUAL "GNU" AND
   CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64" AND
   SCLOUDPLUS_BACKEND STREQUAL "NEON")
  list(APPEND SCLOUDPLUS_EXTRA_COMPILE_OPTIONS -mtune=tsv110)
endif()

set(SCLOUDPLUS_INSTANCE_SUFFIX "-${SCLOUDPLUS_BACKEND}")
set(SCLOUDPLUS_CANONICAL_INSTANCE "Scloudplus-${SCLOUDPLUS_LEVEL}-${SCLOUDPLUS_FAMILY}-packed10")
set(SCLOUDPLUS_SELECTED_INSTANCE "${SCLOUDPLUS_CANONICAL_INSTANCE}${SCLOUDPLUS_INSTANCE_SUFFIX}")
message(STATUS "Configuring ${SCLOUDPLUS_TIER}/${SCLOUDPLUS_SELECTED_INSTANCE}")

set(SCLOUDPLUS_KEM_INSTANCE_SOURCE "${CMAKE_CURRENT_LIST_DIR}/KEM_AlgorithmInstance.c")
set(SCLOUDPLUS_TEST_SOURCE "${SCLOUDPLUS_TEST_DIR}/test_scloudplus.c")
set(SCLOUDPLUS_KEM_LOOP_SOURCE "${SCLOUDPLUS_TEST_DIR}/test_kem_loop.c")
set(SCLOUDPLUS_TAMPER_SOURCE "${SCLOUDPLUS_TEST_DIR}/test_tamper.c")
set(SCLOUDPLUS_BENCHMARK_SOURCE "${SCLOUDPLUS_BENCH_DIR}/kem_benchmark.c")
set(SCLOUDPLUS_KAT_RANDOM_SOURCE "${SCLOUDPLUS_API_PKC_DIR}/kat_random.c")
set(SCLOUDPLUS_KAT_KEM_SOURCE "${SCLOUDPLUS_API_PKC_DIR}/KAT_KEM.c")
set(SCLOUDPLUS_VERIFY_KAT_SOURCE "${SCLOUDPLUS_TEST_DIR}/verify_kat_kem.c")
set(SCLOUDPLUS_DRNG_SOURCE "${SCLOUDPLUS_API_PKC_DIR}/drng.c")

add_executable(test_scloudplus
  ${SCLOUDPLUS_CORE_SOURCES}
  ${SCLOUDPLUS_KEM_INSTANCE_SOURCE}
  ${SCLOUDPLUS_TEST_SOURCE}
)
add_executable(kem_loop_scloudplus
  ${SCLOUDPLUS_CORE_SOURCES}
  ${SCLOUDPLUS_KEM_INSTANCE_SOURCE}
  ${SCLOUDPLUS_KEM_LOOP_SOURCE}
)
add_executable(tamper_scloudplus
  ${SCLOUDPLUS_CORE_SOURCES}
  ${SCLOUDPLUS_KEM_INSTANCE_SOURCE}
  ${SCLOUDPLUS_TAMPER_SOURCE}
)
add_executable(bench_scloudplus
  ${SCLOUDPLUS_CORE_SOURCES}
  ${SCLOUDPLUS_KEM_INSTANCE_SOURCE}
  ${SCLOUDPLUS_BENCHMARK_SOURCE}
)
add_library(scloudplus_kem SHARED
  ${SCLOUDPLUS_CORE_SOURCES}
  ${SCLOUDPLUS_KEM_INSTANCE_SOURCE}
)
set_target_properties(scloudplus_kem PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  OUTPUT_NAME "scloudplus_kem"
)
add_executable(kat_kem_scloudplus
  ${SCLOUDPLUS_KAT_CORE_SOURCES}
  ${SCLOUDPLUS_KAT_RANDOM_SOURCE}
  ${SCLOUDPLUS_KEM_INSTANCE_SOURCE}
  ${SCLOUDPLUS_KAT_KEM_SOURCE}
  ${SCLOUDPLUS_DRNG_SOURCE}
)
add_executable(verify_kat_kem
  ${SCLOUDPLUS_KAT_CORE_SOURCES}
  ${SCLOUDPLUS_KAT_RANDOM_SOURCE}
  ${SCLOUDPLUS_KEM_INSTANCE_SOURCE}
  ${SCLOUDPLUS_VERIFY_KAT_SOURCE}
  ${SCLOUDPLUS_DRNG_SOURCE}
)

foreach(target
    test_scloudplus
    kem_loop_scloudplus
    tamper_scloudplus
    bench_scloudplus
    kat_kem_scloudplus
    verify_kat_kem
    scloudplus_kem)
  target_include_directories(${target} PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}
    ${SCLOUDPLUS_CORE_DIR}/include
    ${SCLOUDPLUS_CORE_DIR}/common
    ${SCLOUDPLUS_API_PKC_DIR}
  )
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -std=c99 -Wall -Wextra -Wpedantic)
    if(SCLOUDPLUS_EXTRA_COMPILE_OPTIONS)
      target_compile_options(${target} PRIVATE ${SCLOUDPLUS_EXTRA_COMPILE_OPTIONS})
    endif()
  endif()
  target_compile_definitions(${target} PRIVATE ${SCLOUDPLUS_EXTRA_COMPILE_DEFINITIONS})
endforeach()
target_compile_definitions(verify_kat_kem PRIVATE SCLOUDPLUS_KAT_INSTANCE_NAME=\"${SCLOUDPLUS_CANONICAL_INSTANCE}\")
