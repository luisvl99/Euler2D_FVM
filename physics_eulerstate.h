#pragma once

// ---------------------------------------------------------------------------
// EulerState
//
//  Stores the four conserved variables of the 2-D Euler equations:
//
//    rho   – density                          [kg/m³]
//    rho_u – x-momentum density               [kg/(m²·s)]
//    rho_v – y-momentum density               [kg/(m²·s)]
//    rho_E – total energy density             [J/m³]
//
//  Arithmetic operators are provided so that flux routines and time-
//  integration loops can be written in natural mathematical notation
//  without temporary helper variables.
// ---------------------------------------------------------------------------

struct EulerState
{
    double rho;
    double rho_u;
    double rho_v;
    double rho_E;


    EulerState()
        : rho(0.0), rho_u(0.0), rho_v(0.0), rho_E(0.0)
    {}


    EulerState(double rho_, double rho_u_, double rho_v_, double rho_E_)
        : rho(rho_), rho_u(rho_u_), rho_v(rho_v_), rho_E(rho_E_)
    {}



    EulerState operator+(const EulerState& o) const
    {
        return { rho   + o.rho,
                rho_u + o.rho_u,
                rho_v + o.rho_v,
                rho_E + o.rho_E };
    }

    EulerState operator-(const EulerState& o) const
    {
        return { rho   - o.rho,
                rho_u - o.rho_u,
                rho_v - o.rho_v,
                rho_E - o.rho_E };
    }

    // scalar * state  and  state * scalar
    EulerState operator*(double s) const
    {
        return { rho*s, rho_u*s, rho_v*s, rho_E*s };
    }

    friend EulerState operator*(double s, const EulerState& U)
    {
        return U * s;
    }

    EulerState operator/(double s) const
    {
        return { rho/s, rho_u/s, rho_v/s, rho_E/s };
    }

    EulerState& operator+=(const EulerState& o)
    {
        rho   += o.rho;
        rho_u += o.rho_u;
        rho_v += o.rho_v;
        rho_E += o.rho_E;
        return *this;
    }

    EulerState& operator-=(const EulerState& o)
    {
        rho   -= o.rho;
        rho_u -= o.rho_u;
        rho_v -= o.rho_v;
        rho_E -= o.rho_E;
        return *this;
    }

    EulerState& operator*=(double s)
    {
        rho   *= s;
        rho_u *= s;
        rho_v *= s;
        rho_E *= s;
        return *this;
    }
};
