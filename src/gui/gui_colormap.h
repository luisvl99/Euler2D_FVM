#pragma once
#include <QColor>
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------------
// Colormap
//
//  Two palettes:
//    Sequential  – jet-like blue→cyan→green→yellow→red  (density, p, Mach…)
//    Diverging   – blue→white→red  (signed fields: u, v)
//
//  Usage:
//    double t = ColorMap::normalize(value, minVal, maxVal);
//    QColor c = ColorMap::sequential(t);
// ---------------------------------------------------------------------------

class ColorMap
{
public:

    // Map value into [0,1] clamped
    static double normalize(double v, double vmin, double vmax)
    {
        if(std::abs(vmax - vmin) < 1e-30)
            return 0.5;
        double t = (v - vmin) / (vmax - vmin);
        return std::clamp(t, 0.0, 1.0);
    }

    // Sequential: jet-style (blue → cyan → green → yellow → red)
    static QColor sequential(double t)
    {
        // Four-segment linear interpolation through 5 anchors
        struct RGB { double r,g,b; };
        static constexpr RGB anchors[5] = {
            {0.0, 0.0, 0.6},   // t=0.00  deep blue
            {0.0, 0.8, 1.0},   // t=0.25  cyan
            {0.0, 0.9, 0.0},   // t=0.50  green
            {1.0, 1.0, 0.0},   // t=0.75  yellow
            {1.0, 0.0, 0.0},   // t=1.00  red
        };
        t = std::clamp(t, 0.0, 1.0);
        double idx = t * 4.0;
        int    lo  = (int)idx;
        if(lo >= 4) lo = 3;
        double f = idx - lo;
        double r = anchors[lo].r + f*(anchors[lo+1].r - anchors[lo].r);
        double g = anchors[lo].g + f*(anchors[lo+1].g - anchors[lo].g);
        double b = anchors[lo].b + f*(anchors[lo+1].b - anchors[lo].b);
        return QColor::fromRgbF(r, g, b);
    }

    // Diverging: blue → white → red  (good for signed fields centered at 0)
    static QColor diverging(double t)
    {
        t = std::clamp(t, 0.0, 1.0);
        double r, g, b;
        if(t < 0.5)
        {
            double f = t * 2.0;           // 0→1 as t goes 0→0.5
            r = f;
            g = f;
            b = 1.0;
        }
        else
        {
            double f = (t - 0.5) * 2.0;  // 0→1 as t goes 0.5→1
            r = 1.0;
            g = 1.0 - f;
            b = 1.0 - f;
        }
        return QColor::fromRgbF(r, g, b);
    }
};
