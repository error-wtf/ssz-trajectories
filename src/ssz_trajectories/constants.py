"""Physical and mathematical constants for the SSZ framework.

Authors: Carmen N. Wrede, Lino P. Casu
"""
import math

PHI = (1.0 + math.sqrt(5.0)) / 2.0
"""Golden ratio phi = (1+sqrt(5))/2."""

LN_PHI = math.log(PHI)
"""ln(phi)."""

c_light = 299_792_458.0
"""Speed of light in m/s."""

RS_DEFAULT = 1.0
"""Default Schwarzschild radius (normalised units)."""

BLEND_LO = 80.0
"""Lower bound of C2 Hermite blend zone (r/r_s)."""

BLEND_HI = 120.0
"""Upper bound of C2 Hermite blend zone (r/r_s)."""
