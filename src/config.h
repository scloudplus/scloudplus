#ifndef _SCLOUDPLUS_CONFIG_H_
#define _SCLOUDPLUS_CONFIG_H_

// Define whether to use OpenSSL
#define USE_OPENSSL 1 // 1 to use OpenSSL, 0 to not use OpenSSL

// Define the operating system type
#define OS_TYPE_LINUX 1
#define OS_TYPE_WINDOWS 2
#define OS_TYPE_MACOS 3

// Automatically detect the operating system
#if defined(_WIN32) || defined(_WIN64)
#define OS_TARGET OS_TYPE_WINDOWS
#elif defined(__APPLE__) || defined(__MACH__)
#define OS_TARGET OS_TYPE_MACOS
#elif defined(__linux__)
#define OS_TARGET OS_TYPE_LINUX
#else
#error "Unknown operating system"
#endif

// Define the processor architecture
#define ARCH_X86 1
#define ARCH_X64 2
#define ARCH_ARM 3
#define ARCH_ARM64 4

// Automatically detect the processor architecture
#if defined(_M_X64) || defined(__x86_64__) || defined(__amd64__)
#define TARGET ARCH_X64
#elif defined(_M_IX86) || defined(__i386__)
#define TARGET ARCH_X86
#elif defined(_M_ARM) || defined(__arm__)
#define TARGET ARCH_ARM
#elif defined(_M_ARM64) || defined(__aarch64__)
#define TARGET ARCH_ARM64
#else
#error "Unknown processor architecture"
#endif

// Define the compiler type
#define COMPILER_GCC 1
#define COMPILER_CLANG 2
#define COMPILER_MSVC 3

// Automatically detect the compiler
#if defined(__clang__)
#define COMPILER COMPILER_CLANG
#elif defined(__GNUC__) || defined(__GNUG__)
#define COMPILER COMPILER_GCC
#elif defined(_MSC_VER)
#define COMPILER COMPILER_MSVC
#else
#error "Unknown compiler"
#endif

#if (OS_TARGET == OS_TYPE_WINDOWS)
#define ALIGN_HEADER(N) __declspec(align(N))
#define ALIGN_FOOTER(N)
#else
#define ALIGN_HEADER(N)
#define ALIGN_FOOTER(N) __attribute__((aligned(N)))
#endif

#endif // CONFIG_H