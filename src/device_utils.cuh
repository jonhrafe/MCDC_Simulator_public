// device_utils.cuh
#pragma once

// // Device-side helper functions as static inline to prevent multiple definitions
// __device__ static inline double3 make_double3(double x, double y, double z) {
//     double3 v;
//     v.x = x;
//     v.y = y;
//     v.z = z;
//     return v;
// }

__device__ static inline double3 operator-(const double3 &a, const double3 &b) {
    return make_double3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__device__ static inline double dot3(const double3 &a, const double3 &b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

__device__ static inline double3 cross3(const double3 &a, const double3 &b) {
    return make_double3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}
