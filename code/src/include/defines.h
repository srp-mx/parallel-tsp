#pragma once

#include <stdint.h>
#include <assert.h>

#define internal static
#define global_variable static
#define local_persist static

#define IGNORE_RESULT(expr) do { if (expr) {} } while (0)

typedef uint8_t u8; static_assert(sizeof(u8) == 1, "u8 not 1 byte");
typedef uint16_t u16; static_assert(sizeof(u16) == 2, "u16 not 2 bytes");
typedef uint32_t u32; static_assert(sizeof(u32) == 4, "u32 not 4 bytes");
typedef uint64_t u64; static_assert(sizeof(u64) == 8, "u64 not 8 bytes");

typedef int8_t i8; static_assert(sizeof(i8) == 1, "i8 not 1 byte");
typedef int16_t i16; static_assert(sizeof(i16) == 2, "i16 not 2 bytes");
typedef int32_t i32; static_assert(sizeof(i32) == 4, "i32 not 4 bytes");
typedef int64_t i64; static_assert(sizeof(i64) == 8, "i64 not 8 bytes");

typedef int32_t b32; static_assert(sizeof(b32) == 4, "b32 not 4 bytes");

typedef float r32; static_assert(sizeof(r32) == 4, "r32 not 4 bytes");
typedef double r64; static_assert(sizeof(r64) == 8, "r64 not 8 bytes");
