/********************************************************************************************
 * ds_benchmark.h: Cross-platform benchmarking macros for C
 *
 * Features:
 *   - Supports Windows, Linux and macOS
 *   - Uses high-resolution timers
 *   - On x86: Provides CPU cycle counting via RDTSC
 *   - Simple statistics: iterations, total time, avg time, cycles
 *   - Two testing modes: fixed iterations or fixed duration
 *
 * Usage:
 *   TIME_OPERATION_ITERATIONS(op, op_name, it) - Run operation 'it' times
 *   TIME_OPERATION_SECONDS(op, op_name, secs)  - Run operation for 'secs' seconds
 *
 * Output:
 *   Operation, Iterations, Total Time (s), Avg Time (us), Total Cycles, Avg Cycles
 *
 * This is free and unencumbered software released into the public domain.
 ********************************************************************************************/

#ifndef _DS_BENCHMARK_H
#define _DS_BENCHMARK_H

#include <stdio.h>
#include <stdint.h>

// Platform detection
#if defined(_WIN32) || defined(_WIN64)
#define OS_WINDOWS

#include <windows.h>

#elif defined(__APPLE__)
#define OS_MACOS
#include <mach/mach_time.h>
#elif defined(__linux__) || defined(__unix__)
#define OS_LINUX
#include <time.h>
#else
#error "Unsupported platform"
#endif

// Architecture detection
#if defined(__x86_64__) || defined(__i386__) || defined(_M_IX86) || defined(_M_X64)
#define ARCH_X86
#endif

// High-resolution timer functions
#if defined(OS_WINDOWS)

static inline uint64_t get_timer_freq() {
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    return freq.QuadPart;
}

static inline uint64_t get_timer_value() {
    LARGE_INTEGER counter;
    QueryPerformanceCounter(&counter);
    return counter.QuadPart;
}

static double timer_to_seconds(uint64_t start, uint64_t end, uint64_t freq) {
    return (double) (end - start) / (double) freq;
}

#elif defined(OS_MACOS)
static inline uint64_t get_timer_freq() {
        mach_timebase_info_data_t timebase;
        mach_timebase_info(&timebase);
        return (uint64_t)(1e9 * (double)timebase.denom / (double)timebase.numer);
    }

    static inline uint64_t get_timer_value() {
        return mach_absolute_time();
    }

    static double timer_to_seconds(uint64_t start, uint64_t end, uint64_t freq) {
        return (double)(end - start) / (double)freq;
    }
#elif defined(OS_LINUX)
    static inline uint64_t get_timer_freq() {
        return 1000000000; // nanoseconds
    }

    static inline uint64_t get_timer_value() {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (uint64_t)ts.tv_sec * 1000000000 + (uint64_t)ts.tv_nsec;
    }

    static double timer_to_seconds(uint64_t start, uint64_t end, uint64_t freq) {
        return (double)(end - start) / (double)freq;
    }
#endif

// CPU cycle counting for x86
#if defined(ARCH_X86)
#if defined(_MSC_VER)
#include <intrin.h>
#define rdtsc() __rdtsc()
#else

static inline uint64_t rdtsc() {
    unsigned int lo, hi;
    __asm__ volatile ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t) hi << 32) | lo;
}

#endif
#endif

// Benchmark macros
#define DEFINE_TIMER_VARIABLES \
    uint64_t _bench_timer_freq = get_timer_freq(); \
    uint64_t _bench_start_time, _bench_end_time; \
    uint64_t _bench_total_time = 0; \
    uint64_t _bench_iterations = 0; \
    uint64_t _bench_total_cycles = 0;

#if defined(ARCH_X86)
#define DEFINE_CYCLE_VARIABLES \
        uint64_t _bench_start_cycles, _bench_end_cycles;

#define START_CYCLE_COUNT \
        _bench_start_cycles = rdtsc();

#define STOP_CYCLE_COUNT \
        _bench_end_cycles = rdtsc(); \
        if (_bench_end_cycles < _bench_start_cycles) { \
            _bench_end_cycles += (uint64_t)1 << 32; \
        } \
        _bench_total_cycles += (_bench_end_cycles - _bench_start_cycles);
#else
#define DEFINE_CYCLE_VARIABLES
#define START_CYCLE_COUNT
#define STOP_CYCLE_COUNT
#endif

#define INITIALIZE_TIMER \
    _bench_total_time = 0; \
    _bench_iterations = 0; \
    _bench_total_cycles = 0;

#define START_TIMER \
    _bench_start_time = get_timer_value(); \
    START_CYCLE_COUNT

#define STOP_TIMER \
    STOP_CYCLE_COUNT \
    _bench_end_time = get_timer_value(); \
    _bench_total_time += (_bench_end_time - _bench_start_time); \
    _bench_iterations++;

#define PRINT_TIMER_HEADER \
    printf("%-30s %12s %15s %15s %15s %15s\n", \
           "Operation", "Iterations", "Total Time (s)", "Avg Time (us)", \
           "Total Cycles", "Avg Cycles");

#define PRINT_TIMER_AVG(op_name) \
    do { \
        double total_seconds = timer_to_seconds(0, _bench_total_time, _bench_timer_freq); \
        double avg_us = (total_seconds * 1e6) / _bench_iterations; \
        double avg_cycles = _bench_iterations ? (double)_bench_total_cycles / _bench_iterations : 0; \
        printf("%-30s %12llu %15.6f %15.3f %15llu %15.0f\n", \
               (op_name), _bench_iterations, total_seconds, avg_us, \
               _bench_total_cycles, avg_cycles); \
    } while (0)

// Fixed iteration count testing
#define TIME_OPERATION_ITERATIONS(op, op_name, it) \
    do { \
        DEFINE_TIMER_VARIABLES \
        DEFINE_CYCLE_VARIABLES \
        INITIALIZE_TIMER \
        for (uint64_t _i = 0; _i < (it); _i++) { \
            START_TIMER; \
            { op; } \
            STOP_TIMER; \
        } \
        PRINT_TIMER_AVG(op_name); \
    } while (0)

// Fixed duration testing
#define TIME_OPERATION_SECONDS(op, op_name, secs) \
    do { \
        DEFINE_TIMER_VARIABLES \
        DEFINE_CYCLE_VARIABLES \
        INITIALIZE_TIMER \
        const double _target_sec = (secs); \
        while (timer_to_seconds(0, _bench_total_time, _bench_timer_freq) < _target_sec) { \
            START_TIMER; \
            { op; } \
            STOP_TIMER; \
        } \
        PRINT_TIMER_AVG(op_name); \
    } while (0)

#endif // _DS_BENCHMARK_H