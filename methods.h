#ifndef CM_METHODS_H
#define CM_METHODS_H
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <math.h>
#include <functional>

using ld = long double;
const ld EPS = 1e-9;

ld simple_iteration_method (const std::function<ld(ld)>& phi_f, const std::function<ld(ld)>& d_phi_f, ld a, ld b, ld epsilon) {
    ld q;
    ld prev = (a + b) * 0.5;
    ld x = prev;

    do {
        prev = x;
        x = phi_f(prev);

        q = abs(d_phi_f(x));
        if (q >= 1) {
            throw std::runtime_error("The convergence condition is not satisfied!");
        }
//        std::cout << q << std::endl;

    } while (((q / (1 - q)) * abs(x - prev)) > epsilon);

    return x;
}

ld dichotomy_method (const std::function<ld(ld)>& f, ld a, ld b, ld epsilon) {
    if (f(a) * f(b) > 0) {
        throw std::runtime_error("The same sign at the ends of segment!");
    }

    ld mid;
    while ((b - a) > 2 * epsilon) {
        mid = (a + b) * 0.5;

        if (abs(f(mid)) < epsilon) {
            return mid;
        }

        if (f(a) * f(mid) < 0) {
            b = mid;
        } else {
            a = mid;
        }
    }

    return (a + b) * 0.5;
}

ld newton_method (const std::function<ld(ld)>& func, const std::function<ld(ld)>& dfunc,
                  const std::function<ld(ld)>& ddfunc, ld a, ld b,ld epsilon) {
    ld x, prev;

    if (func(a) * func(b) > 0) {
        throw std::runtime_error("The same sign at the ends of segment!");
    }

    if (func(a) * ddfunc(a) > 0) {
        x = a;
    } else if (func(b) * ddfunc(b) > 0) {
        x = b;
    } else {
        x = (a + b) * 0.5;
    }

    ld df;
    do {
        prev = x;

        df = dfunc(x);

        if (abs(df) < epsilon) {
            throw std::runtime_error("Division by zero!");
        }

        x = x - func(x) / df;
    } while (std::abs (prev - x) > epsilon);

    return x;
}

ld secant_method (const std::function<ld(ld)>& func, const std::function<ld(ld)>& ddfunc,
                  ld epsilon, ld a, ld b) {
    ld x, prev, next;

    if (func(a) * func(b) > 0) {
        throw std::runtime_error("The same sign at the ends of segment!");
    }

    if (func(a) * ddfunc(a) > 0) {
        prev = a;
    } else if (func(b) * ddfunc(b) > 0) {
        prev = b;
    } else {
        prev = (a + b) * 0.5;
    }

    x = prev + epsilon;

    ld fx, fxp;
    do {
        fx = func(x);
        fxp = func(prev);

        next = x - (fx * (x - prev)) / (fx - fxp);

        prev = x;
        x = next;
    } while (abs(prev - x) > epsilon);

    return x;
}

#endif //CM_METHODS_H
