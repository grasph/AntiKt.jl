# Copyright (c) 2022, Philippe Gras
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

"""Provides anti-kt hadronic jet clustering based on a Fastjet algorithm
"""
module AntiKt

export PseudoJet, antikt

include("pseudojet.jl")

include("antiktalgo.jl")

antikt = _faster_tiled_N2_cluster

end
