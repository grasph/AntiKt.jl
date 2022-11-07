# Copyright (c) 2022, Philippe Gras
#
# Adapted from PseudoJet class of c++ code from the Fastjet
# software  (https://fastjet.fr,  hep-ph/0512210,  arXiv:1111.6097)
#
#   Copyright (c) 2005-2020, Matteo Cacciari, Gavin P. Salam and Gregory
#   Soyez
#
#
#----------------------------------------------------------------------
# This file is part of AntiKt.jl.
#
#  AntiKt.jl is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  The algorithms that underlie FastJet have required considerable
#  development. They are described in the original FastJet paper,
#  hep-ph/0512210 and in the manual, arXiv:1111.6097. If you use
#  FastJet as part of work towards a scientific publication, please
#  quote the version you use and include a citation to the manual and
#  optionally also to hep-ph/0512210.
#
#  AntiKet.jl is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with FastJet. If not, see <http:#www.gnu.org/licenses/>.
#----------------------------------------------------------------------


"""Interface for composite types that includes fields px, py, py, and E
that represents the components of a four-momentum vector."""
abstract type FourMomentum end

# Used to protect against parton-level events where pt can be zero
# for some partons, giving rapidity=infinity. KtJet fails in those cases.
const _MaxRap = 1e5

# @ingroup basic_classes
# \class PseudoJet
# Class to contain pseudojets, including minimal information of use to
# jet-clustering routines.
mutable struct PseudoJet<:FourMomentum
    # construct a pseudojet from explicit components
    px::Float64
    py::Float64
    pz::Float64
    E::Float64
    _cluster_hist_index::Int
    _rap::Float64
    _phi::Float64
    _pt2::Float64
    _inv_pt2::Float64
end

PseudoJet(px, py, pz, E) = begin
    pt2 = px*px +py*py
    inv_pt2 = inv(pt2)
    PseudoJet(px, py, pz, E, 0, NaN, NaN, pt2, inv_pt2)
end


import Base.getproperty
Base.getproperty(x::PseudoJet, sym::Symbol) = begin
    if sym === :pt2
        return x._pt2
    elseif sym === :pt
        return sqrt(x._pt2)
    elseif sym == :rap
        return _get_rap(x)
    elseif sym == :phi
        return _get_phi(x)
    else
        return getfield(x, sym)
    end
end


set_momentum(j::PseudoJet, px, py, pz, E) = begin
    j.px = px
    j.py = py
    j.pz = pz
    j.E = E
    j._pt2 = px^2 + py^2
    j._inv_pt2 = 1.0/j._pt2
    j._rap = NaN
    j._phi = NaN
end

_fix_phi(phi::Real) = begin
    if phi < 0.0
        phi + 2π
    elseif phi >= 2π
       phi - 2π  # can happen if phi=-|eps<1e-15|?
    else
        phi
   end
end

_get_phi(p::PseudoJet) = begin
    _fix_phi(p._pt2 == 0.0 ? 0.0 : atan(p.py,p.px))
end

_get_rap(p::PseudoJet) = begin
    if p.E == abs(p.pz) && iszero(p._pt2)
        # Point has infinite rapidity -- convert that into a very large
        #    number, but in such a way that different 0-pt momenta will have
        #    different rapidities (so as to lift the degeneracy between
        #                         them) [this can be relevant at parton-level]
        MaxRapHere = _MaxRap + abs(p.pz)
        p.pz >= 0.0 ?  MaxRapHere : -MaxRapHere
    else
        # get the rapidity in a way that's modestly insensitive to roundoff
        # error when things pz,E are large (actually the best we can do without
        # explicit knowledge of mass)
        effective_m2 = max(0.0, m2(p)) # force non tachyonic mass
        E_plus_pz    = p.E + abs(p.pz) # the safer of p+, p-
        rap_tmp = 0.5*log((p._pt2 + effective_m2)/(E_plus_pz*E_plus_pz))
        p.pz > 0 ? - rap_tmp : rap_tmp
    end
end

phi(p::PseudoJet) = phi_02pi(p)

phi_02pi(p::PseudoJet) = p.phi

rap(p::PseudoJet) = p.rap

pt2(p::PseudoJet) = p.pt2

eta(p::PseudoJet) = p.rap

"Returns the scalar transverse momentum"
pt(p::PseudoJet) = sqrt(p._pt2)

"Returns the squared invariant mass"
m2(p::PseudoJet) = (p.E + p.pz)*(p.E-p.pz) - p._pt2

# returns the invariant mass
# (If m2() is negative then -sqrt(-m2()) is returned, as in CLHEP)
m(p::PseudoJet) = begin
    _m = sqrt(p.m2())
    _m < 0. ? -sqrt(-m) : sqrt(m)
end

import Base.+

+(j1::PseudoJet, j2::PseudoJet) = PseudoJet(j1.px + j2.px, j1.py + j2.py,
                                            j1.pz + j2.pz, j1.E + j2.E)
