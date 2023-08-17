# Copyright (c) 2022, Philippe Gras
#
# Adapted from ClusterSequence_Tiled_N2.cc c++ code from the Fastjet
# software  (https://fastjet.fr,  hep-ph/0512210,  arXiv:1111.6097)
#
#   Copyright (c) 2005-2020, Matteo Cacciari, Gavin P. Salam and Gregory
#   Soyez
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
const Invalid=-3
const InexistentParent=-2
const BeamJet=-1

import Base.convert

mydebug(x...) = println(">>+> ", x...)

"""for debugging to compare with indices used in the C++ implementation
"""
flatindex(c::CartesianIndex{2},a::Matrix) = (c[1]-1)*size(a)[2] + c[2]

struct HistoryElement
    """Index in _history where first parent of this jet
    was created (InexistentParent if this jet is an
    original particle)"""
    parent1::Int

    """index in _history where second parent of this jet
    was created (InexistentParent if this jet is an
    original particle); BeamJet if this history entry
    just labels the fact that the jet has recombined
    with the beam)"""
    parent2::Int

    """index in _history where the current jet is
    recombined with another jet to form its child. It
    is Invalid if this jet does not further
    recombine."""
    child::Int

    """index in the _jets vector where we will find the
    PseudoJet object corresponding to this jet
    (i.e. the jet created at this entry of the
    history). NB: if this element of the history
    corresponds to a beam recombination, then
    jetp_index=Invalid."""
    jetp_index::Int

    """the distance corresponding to the recombination
       at this stage of the clustering."""
    dij::Float64

    """the largest recombination distance seen
       so far in the clustering history."""
    max_dij_so_far::Float64
end

HistoryElement(jetp_index) = HistoryElement(InexistentParent, InexistentParent, Invalid, jetp_index, 0.0, 0.0)


"""Structure analogous to BriefJet, but with the extra information
needed for dealing with tiles"""
mutable struct TiledJet
    id::Int
    eta::Float64
    phi::Float64
    kt2::Float64
    NN_dist::Float64

    jets_index::Int
    tile_index::Int
    diJ_posn::Int


    "Nearest neighbour"
    NN::TiledJet

    previous::TiledJet
    next::TiledJet
    TiledJet(::Type{Nothing}) = begin
        t = new(-1, 0., 0., 0., 0., -1, 0, 0)
        t.NN = t.previous = t.next = t
        t
    end
    
    TiledJet(id, eta, phi, kt2, NN_dist,
             jet_index, tile_index, diJ_posn,
             NN, previous, next) = new(id, eta, phi, kt2, NN_dist,
                                       jet_index, tile_index, diJ_posn,
                                       NN, previous, next)
end


const noTiledJet::TiledJet = TiledJet(Nothing)

isvalid(t::TiledJet) = !(t === noTiledJet)

""" Move a TiledJet in front of a TiledJet list element.

The jet to move can be an isolated jet, a jet from another list or a jet from the same list.

"""
insert!(nextjet::TiledJet, jettomove::TiledJet) = begin
    if !isnothing(nextjet)
        nextjet.previous  = jettomove
    end

    jettomove.next = nextjet
    jettomove.previous = nextjet.previous
    nextjet = jettomove
end

""" Detach a TiledJet from its list.
"""
detach!(jet::TiledJet) = begin
    if !isnothing(jet.previous)
        jet.previous.next = jet.next
    end
    if !isnothing(jet.next)
        jet.next.previous = jet.previous
    end
    jet.next = jet.previous = noTiledJet
end


TiledJet(id) = TiledJet(id, 0., 0., 0., 0.,
                        0, 0, 0,
                        noTiledJet, noTiledJet, noTiledJet)

import Base.copy
copy(j::TiledJet) = TiledJet(j.id, j.eta, j.phi, j.kt2, j.NN_dist, j.jets_index, j.tile_index, j.diJ_posn, j.NN, j.previous, j.next)

"""Computes distance in the (eta,phi)-plane
between two jets."""
_tj_dist(jetA, jetB) = begin
    dphi = π - abs(π - abs(jetA.phi - jetB.phi))
    deta = jetA.eta - jetB.eta
    return dphi*dphi + deta*deta
end

_tj_diJ(jet) = begin
    kt2 = jet.kt2
    if isvalid(jet.NN) && jet.NN.kt2 < kt2
        kt2 = jet.NN.kt2
    end
    return jet.NN_dist * kt2
end


const _n_tile_center = 1
const _n_tile_left_neighbours = 4
const _tile_right_neigbour_indices = 6:9
const _n_tile_right_neighbours = 4
const _n_tile_neighbours = 9

const neigh_init = fill(nothing,_n_tile_neighbours)

struct TilingDef{F,I}
    _tiles_eta_min::F
    _tiles_eta_max::F
    _tile_size_eta::F
    _tile_size_phi::F
    _n_tiles_eta::I
    _n_tiles_phi::I
    _n_tiles::I
    _tiles_ieta_min::I
    _tiles_ieta_max::I
end

const tile_left = -1
const tile_central = 0
const tile_right = 1

struct Tiling{F,I}
    setup::TilingDef{F,I}
    tiles::Matrix{TiledJet}
    positions::Matrix{Int}
    tags::Matrix{Bool}
end

Tiling(setup::TilingDef) = begin
    t = Tiling(setup,
               fill(noTiledJet, (setup._n_tiles_eta, setup._n_tiles_phi)),
               fill(0, (setup._n_tiles_eta, setup._n_tiles_phi)),
               fill(false, (setup._n_tiles_eta, setup._n_tiles_phi)))
    @inbounds for iphi = 1:setup._n_tiles_phi
        #order of the following two statement is important
        #to have position = tile_right in case n_tiles_eta = 1
        t.positions[1, iphi] = tile_left
        t.positions[setup._n_tiles_eta, iphi] = tile_right
    end
    t
end

struct Surrounding{N}
    indices::NTuple{N, Int}
end

import Base.iterate

Base.iterate(x::T) where {T<:Surrounding} = (x.indices[1], 2)
Base.iterate(x::Surrounding{0}) = nothing
Base.iterate(x::Surrounding{1}, state) = nothing
Base.iterate(x::Surrounding{2}, state) = nothing
Base.iterate(x::Surrounding{3}, state) = state > 3 ? nothing : (x.indices[state], state+1)
Base.iterate(x::Surrounding{4}, state) = state > 4 ? nothing : (x.indices[state], state+1)
Base.iterate(x::Surrounding{6}, state) = state > 6 ? nothing : (x.indices[state], state+1)
Base.iterate(x::Surrounding{9}, state) = state > 9 ? nothing : (x.indices[state], state+1)

import Base.length
length(x::Surrounding{N}) where N = N


surrounding(center::Int, tiling::Tiling) = begin
    #                        4|6|9
    #                        3|1|8
    #                        2|5|7
    #  -> η

    iphip = mod1(center + tiling.setup._n_tiles_eta, tiling.setup._n_tiles)
    iphim = mod1(center - tiling.setup._n_tiles_eta, tiling.setup._n_tiles)

    if tiling.setup._n_tiles_eta == 1
        return Surrounding{3}((center, iphim, iphip))
    elseif tiling.positions[center] == tile_right
        return Surrounding{6}((center, iphim, iphip, iphim - 1, center - 1, iphip - 1))
    elseif tiling.positions[center] == tile_central
        return Surrounding{9}((center, iphim - 1, center - 1, iphip - 1,
                               iphim, iphip,
                               iphim + 1, center + 1, iphip + 1))
    else #tile_left
        return Surrounding{6}((center, iphim, iphip,
                               iphim + 1, center + 1, iphip + 1))
    end
end

rightneighbours(center::Int, tiling::Tiling) = begin
    #                         |1|4
    #                         | |3
    #                         | |2
    #  -> η

    iphip = mod1(center + tiling.setup._n_tiles_eta, tiling.setup._n_tiles)
    iphim = mod1(center - tiling.setup._n_tiles_eta, tiling.setup._n_tiles)

    if tiling.positions[center] == tile_right
        return Surrounding{1}((iphip,))
    else
        return Surrounding{4}((iphip, iphim + 1, center + 1, iphip + 1))
    end
end


struct ClusterSequence
    """This contains the physical PseudoJets; for each PseudoJet one
can find the corresponding position in the _history by looking
at _jets[i].cluster_hist_index()"""
    jets::Vector{PseudoJet}

    """This vector will contain the branching history; for each stage,
_history[i].jetp_index indicates where to look in the _jets
vector to get the physical PseudoJet."""
    history::Vector{HistoryElement}

    """PseudoJet tiling"""
    tiling::Tiling{Float64,Int}

    """Total energy of the event"""
    Qtot
end

##FIXME: Use OffsetArrays
#"""Reasonably robust return of tile index given ieta and iphi, in particular
#it works even if iphi is negative."""
#_tile_index(tiling_setup, ieta::Integer, iphi::Integer) = begin
#    1 + (ieta-tiling_setup._tiles_ieta_min)*tiling_setup._n_tiles_phi + mod(iphi, tiling_setup._n_tiles_phi)
#end


"""Return the tile index corresponding to the given eta,phi point"""
_tile_index(tiling_setup, eta::Float64, phi::Float64) = begin
    if eta <= tiling_setup._tiles_eta_min
        ieta = 1
    elseif eta >= tiling_setup._tiles_eta_max
        ieta = tiling_setup._n_tiles_eta
    else
        ieta = 1 + unsafe_trunc(Int, (eta - tiling_setup._tiles_eta_min) / tiling_setup._tile_size_eta)
        # following needed in case of rare but nasty rounding errors
        if ieta > tiling_setup._n_tiles_eta
            ieta = tiling_setup._n_tiles_eta
        end
    end
    # allow for some extent of being beyond range in calculation of phi
    # as well
    #iphi = (int(floor(phi/_tile_size_phi)) + _n_tiles_phi) % _n_tiles_phi;
    # with just int and no floor, things run faster but beware
    #iphi = mod(unsafe_trunc(Int, (phi + 2π) / tiling_setup._tile_size_phi),
    #           tiling_setup._n_tiles_phi)
    iphi = min(unsafe_trunc(Int, phi  / tiling_setup._tile_size_phi), tiling_setup._n_tiles_phi)
    return iphi * tiling_setup._n_tiles_eta + ieta
end

""""Have a binning of rapidity that goes from -nrap to nrap
in bins of size 1; the left and right-most bins include
include overflows from smaller/larger rapidities"""
determine_rapidity_extent(particles) = begin

    length(particles) == 0 && return 0., 0.

    nrap = 20
    nbins = 2*nrap
    counts = zeros(nbins)

    # get the minimum and maximum rapidities and at the same time bin
    # the multiplicities as a function of rapidity to help decide how
    # far out it's worth going
    minrap =  floatmax(Float64)
    maxrap = -floatmax(Float64)

    ibin = 0
    for p in particles
        # ignore particles with infinite rapidity
        p.E != abs(p.pz) || continue

        y = rap(p)
        minrap = min(minrap, y)
        maxrap = max(maxrap, y)

        # now bin the rapidity to decide how far to go with the tiling.
        # Remember the bins go from ibin=1 (rap=-infinity..-19)
        # to ibin = nbins (rap=19..infinity for nrap=20)
        ibin = clamp(1 + unsafe_trunc(Int, y + nrap), 1, nbins)
        @inbounds counts[ibin] += 1
    end

    # now figure out the particle count in the busiest bin
    max_in_bin = maximum(counts)

    # and find minrap, maxrap such that edge bin never contains more
    # than some fraction of busiest, and at least a few particles; first do
    # it from left. NB: the thresholds chosen here are largely
    # guesstimates as to what might work.
    #
    # 2014-07-17: in some tests at high multiplicity (100k) and particles going up to
    #             about 7.3, anti-kt R=0.4, we found that 0.25 gave 20% better run times
    #             than the original value of 0.5.
    allowed_max_fraction = 0.25

    # the edge bins should also contain at least min_multiplicity particles
    min_multiplicity = 4

    # now calculate how much we can accumulate into an edge bin
    allowed_max_cumul = floor(max(max_in_bin * allowed_max_fraction,
                                  min_multiplicity))

    # make sure we don't require more particles in a bin than max_in_bin
    allowed_max_cumul = min(max_in_bin, allowed_max_cumul)

    # start scan over rapidity bins from the left, to find out minimum rapidity of tiling
    cumul_lo = 0.0
    cumul2 = 0.0
    ibin_lo = 1
    while ibin_lo <= nbins
        @inbounds cumul_lo += counts[ibin_lo]
        if cumul_lo >= allowed_max_cumul
            minrap = max(minrap, ibin_lo - nrap - 1)
            break
        end
        ibin_lo += 1
    end
    @assert ibin_lo != nbins # internal consistency check that you found a bin
    cumul2 += cumul_lo^2

    # then do it from right, to find out maximum rapidity of tiling
    cumul_hi = 0.0
    ibin_hi = nbins
    while ibin_hi >= 1
        @inbounds cumul_hi += counts[ibin_hi]
        if cumul_hi >= allowed_max_cumul
            maxrap = min(maxrap, ibin_hi - nrap)
            break
        end
        ibin_hi -= 1
    end

    @assert ibin_hi >= 0 # internal consistency check that you found a bin

    # consistency check
    @assert ibin_hi >= ibin_lo

    # now work out cumul2
    if ibin_hi == ibin_lo
        # if there is a single bin (potentially including overflows
        # from both sides), cumul2 is the square of the total contents
        # of that bin, which we obtain from cumul_lo and cumul_hi minus
        # the double counting of part that is contained in both
        # (putting double)
        @inbounds cumul2 = (cumul_lo + cumul_hi - counts[ibin_hi]) ^ 2
    else
        # otherwise we have a straightforward sum of squares of bin
        # contents
        cumul2 += cumul_hi^2
    end

    # now get the rest of the squared bin contents
    for ibin in (ibin_lo + 1):ibin_hi
        @inbounds cumul2 += counts[ibin]^2
    end

    minrap, maxrap
end

"""Remove a jet from a tiling"""
_tj_remove_from_tiles!(tiling, jet) = begin
    if !isvalid(jet.previous)
        # we are at head of the tile, so reset it.
        # If this was the only jet on the tile then tile->head will now be NULL
        tiling.tiles[jet.tile_index] = jet.next
    else
        # adjust link from previous jet in this tile
        jet.previous.next = jet.next
    end

    if isvalid(jet.next)
        # adjust backwards-link from next jet in this tile
        jet.next.previous = jet.previous
    end

    jet.next = jet.previous = noTiledJet #to be clean, but not mandatory
end

#----------------------------------------------------------------------
#/ Set up the tiles:
#/  - decide the range in eta
#/  - allocate the tiles
#/  - set up the cross-referencing info between tiles
#/
#/ The neighbourhood of a tile is set up as follows
#/
#/	        LRR
#/          LXR
#/          LLR
#/
#/ such that tiles is an array containing XLLLLRRRR with pointers
#/                                         |   \ RH_tiles
#/                                         \ surrounding_tiles
#/
#/ with appropriate precautions when close to the edge of the tiled
#/ region.
#/
_initial_tiling(particles, Rparam) = begin

    # first decide tile sizes (with a lower bound to avoid huge memory use with
    # very small R)
    tile_size_eta = max(0.1, Rparam)

    # it makes no sense to go below 3 tiles in phi -- 3 tiles is
    # sufficient to make sure all pair-wise combinations up to pi in
    # phi are possible
    n_tiles_phi   = max(3, floor(Int, 2π/tile_size_eta))

    tile_size_phi = 2π / n_tiles_phi # >= Rparam and fits in 2pi

    tiles_eta_min, tiles_eta_max = determine_rapidity_extent(particles)

    # now adjust the values
    tiles_ieta_min = floor(Int, tiles_eta_min/tile_size_eta)
    tiles_ieta_max = floor(Int, tiles_eta_max/tile_size_eta) #FIXME shouldn't it be ceil ?
    tiles_eta_min = tiles_ieta_min * tile_size_eta
    tiles_eta_max = tiles_ieta_max * tile_size_eta
    n_tiles_eta = tiles_ieta_max - tiles_ieta_min + 1

    tiling_setup = TilingDef(tiles_eta_min, tiles_eta_max,
                             tile_size_eta, tile_size_phi,
                             n_tiles_eta, n_tiles_phi, n_tiles_eta*n_tiles_phi,
                             tiles_ieta_min, tiles_ieta_max)

    # allocate the tiles
    tiling_size = (n_tiles_eta, n_tiles_phi)
    Tiling(tiling_setup)
end


#----------------------------------------------------------------------
# Set up a TiledJet
_tj_set_jetinfo!(jet::TiledJet, cs::ClusterSequence, jets_index, R2) = begin
    @inbounds jet.eta  = rap(cs.jets[jets_index])
    @inbounds jet.phi  = phi_02pi(cs.jets[jets_index])
    @inbounds jet.kt2  = pt2(cs.jets[jets_index]) > 1.e-300 ? 1. / pt2(cs.jets[jets_index]) : 1.e300
    jet.jets_index = jets_index
    # initialise NN info as well
    jet.NN_dist = R2
    jet.NN      = noTiledJet

    # Find out which tile it belonds to
    jet.tile_index = _tile_index(cs.tiling.setup, jet.eta, jet.phi)

    # Insert it into the tile's linked list of jets
    jet.previous = noTiledJet
    @inbounds jet.next = cs.tiling.tiles[jet.tile_index]
    if isvalid(jet.next) jet.next.previous = jet; end
    @inbounds cs.tiling.tiles[jet.tile_index] = jet
    nothing
end


Base.iterate(tj::TiledJet) = (tj, tj)
Base.iterate(tj::TiledJet, state::TiledJet) = begin
    isvalid(state.next) ? (state.next::TiledJet, state.next::TiledJet) : nothing
end

# #----------------------------------------------------------------------
# #/ output the contents of the tiles
# _print_tiles(tiling::Tiling) = begin
#   for (itile, tile) in enumerate(tiling.tiles)
#     print("Tile ", itile, " = ")
#     print(join(map(x->x.jets_index, tile.jets, " ")))
#     prinln()
#   end
# end


#----------------------------------------------------------------------
"""Adds to the vector tile_union the tiles that are in the neighbourhood
of the specified tile_index, including itself and whose tagged status are
false ---start adding from position n_near_tiles-1, and increase n_near_tiles as
you go along. When a neighbour is added its tagged status is set to true.

Returns the updated number of near_tiles."""
function _add_untagged_neighbours_to_tile_union(center_index,
                                                tile_union, n_near_tiles,
                                                tiling)
    for tile_index in surrounding(center_index, tiling)
        @inbounds if !tiling.tags[tile_index]
            n_near_tiles += 1
            tile_union[n_near_tiles] = tile_index
            tiling.tags[tile_index] = true
        else
        end
    end
    n_near_tiles
end

#----------------------------------------------------------------------
# Initialise the clustering history in a standard way,
# Takes as input the list of stable particles as input
# Returns the history and the total event energy.
_initial_history(particles) = begin

    # reserve sufficient space for everything
    history = Vector{HistoryElement}(undef, length(particles))
    sizehint!(history, 2*length(particles)) #FIXME does it bring any significant performance improvement?

    Qtot = 0.

    for i in eachindex(particles)
        history[i] = HistoryElement(i)

        # get cross-referencing right from PseudoJets
        particles[i]._cluster_hist_index = i

        # determine the total energy in the event
        Qtot += particles[i].E
    end
    history, Qtot
end


#----------------------------------------------------------------------
# initialise the history in a standard way
_add_step_to_history!(cs::ClusterSequence, parent1, parent2, jetp_index, dij) = begin
    max_dij_so_far = max(dij, cs.history[end].max_dij_so_far)
    push!(cs.history, HistoryElement(parent1, parent2, #=child=# Invalid,
                                     jetp_index, dij, max_dij_so_far))

    local_step = length(cs.history)

    ##ifndef __NO_ASSERTS__
    #assert(local_step == step_number);
    ##endif

    # sanity check: make sure the particles have not already been recombined
    #
    # Note that good practice would make this an assert (since this is
    # a serious internal issue). However, we decided to throw an
    # InternalError so that the end user can decide to catch it and
    # retry the clustering with a different strategy.

    @assert parent1 >= 1
    if cs.history[parent1].child != Invalid
        throw(ErrorException("Internal error. Trying to recombine an object that has previsously been recombined."))
    end

    hist_elem = cs.history[parent1]
    cs.history[parent1] = @set hist_elem.child = local_step

    if parent2 >= 1
        cs.history[parent2].child == Invalid || error("Internal error. Trying to recombine an object that has previsously been recombined.  Parent " * string(parent2) * "'s child index " * string(cs.history[parent1].child) * ". Parent jet index: " * string(cs.history[parent2].jetp_index) * ".")
        hist_elem = cs.history[parent2]
        cs.history[parent2] = @set hist_elem.child = local_step
    end

    # get cross-referencing right from PseudoJets
    if jetp_index != Invalid
        @assert jetp_index >= 1
        cs.jets[jetp_index]._cluster_hist_index = local_step
    end

    #if (_writeout_combinations) {
    #  cout << local_step << ": "
    #	 << parent1 << " with " << parent2
    #	 << "; y = "<< dij<<endl;
    #}
end


"""Carries out the bookkeeping associated with the step of recombining
jet_i and jet_j (assuming a distance dij) and returns the index
of the recombined jet, newjet_k."""
_do_ij_recombination_step!(cs::ClusterSequence, jet_i, jet_j, dij) = begin

    # Create the new jet by recombining the first two with
    # the E-scheme
    #
    push!(cs.jets, cs.jets[jet_i] + cs.jets[jet_j])

    # get its index
    newjet_k = length(cs.jets)

    # get history index
    newstep_k = length(cs.history) + 1

    # and provide jet with the info
    cs.jets[newjet_k]._cluster_hist_index = newstep_k

    # finally sort out the history
    hist_i = cs.jets[jet_i]._cluster_hist_index
    hist_j = cs.jets[jet_j]._cluster_hist_index

    _add_step_to_history!(cs, minmax(hist_i, hist_j)...,
                         newjet_k, dij)

    newjet_k
end

"""Carries out the bookkeeping associated with the step of recombining
jet_i with the beam"""
_do_iB_recombination_step!(cs::ClusterSequence, jet_i, diB) = begin
    # recombine the jet with the beam
    _add_step_to_history!(cs, cs.jets[jet_i]._cluster_hist_index, BeamJet,
                         Invalid, diB)
end
#----------------------------------------------------------------------
# return all inclusive jets of a ClusterSequence with pt > ptmin
inclusive_jets(cs::ClusterSequence, ptmin = 0.) = begin
    dcut = ptmin*ptmin
    jets_local = PseudoJet[]
    sizehint!(jets_local, length(cs.jets))
    # For inclusive jets with a plugin algorithm, we make no
    # assumptions about anything (relation of dij to momenta,
    # ordering of the dij, etc.)
    for elt in Iterators.reverse(cs.history)
        elt.parent2 == BeamJet || continue
        iparent_jet = cs.history[elt.parent1].jetp_index
        jet = cs.jets[iparent_jet]
        if pt2(jet) >= dcut
            push!(jets_local, jet)
        end
    end
    jets_local
end


struct _DiJ_plus_link
    diJ::Float64  # the distance
    jet::TiledJet # the jet (i) for which we've found this distance
end

find_best(diJ, n) = begin
    best = 1
    @inbounds diJ_min = diJ[1]
    @turbo for here in 2:n
        newmin = diJ[here] < diJ_min
        best = newmin ? here : best
        diJ_min = newmin ? diJ[here] : diJ_min
    end
    diJ_min, best
end


#----------------------------------------------------------------------
#/ run a tiled clustering
_faster_tiled_N2_cluster(particles, Rparam, ptmin = 0.0) = begin
    R2 = Rparam * Rparam
    invR2 = 1.0/R2


    # this will ensure that we can point to jets without difficulties
    # arising.
    jets = PseudoJet[]
    sizehint!(jets, length(particles)*2)
    resize!(jets, length(particles))

    # insert initial jets this way so that any type that can be
    # converted to a pseudojet will work fine
    copyto!(jets, particles)

    history, Qtot = _initial_history(jets)

    tiling = _initial_tiling(particles, Rparam)

    cs = ClusterSequence(jets, history, tiling, Qtot)

    tiledjets = similar(cs.jets, TiledJet)

    # will be used quite deep inside loops, but declare it here so that
    # memory (de)allocation gets done only once
    tile_union = Vector{Int}(undef, 3*_n_tile_neighbours)

    # initialise the basic jet info
    for ij in eachindex(tiledjets)
        tiledjets[ij] = TiledJet(ij)
        _tj_set_jetinfo!(tiledjets[ij], cs, ij, R2)
    end

    #_tj_set_jetinfo! links the elements, such that
    #iterating on tiledjet[0] will pass through all the
    #jets

    # set up the initial nearest neighbour information
    for tile in cs.tiling.tiles
        isvalid(tile) || continue
        for jetA in tile
            for jetB in tile
                if jetB == jetA break; end
                dist = _tj_dist(jetA, jetB)
                if (dist < jetA.NN_dist)
                    jetA.NN_dist = dist
                    jetA.NN = jetB
                end
                if dist < jetB.NN_dist
                    jetB.NN_dist = dist
                    jetB.NN = jetA
                end
            end #next jetB
        end #next jetA

        # look for neighbour jets n the neighbour tiles
        for rtile_index in rightneighbours(tile.tile_index, tiling)
            for jetA in tile
                for jetB in @inbounds tiling.tiles[rtile_index]
                    dist = _tj_dist(jetA, jetB)
                    if (dist < jetA.NN_dist)
                        jetA.NN_dist = dist
                        jetA.NN = jetB
                    end
                    if dist < jetB.NN_dist
                        jetB.NN_dist = dist
                        jetB.NN = jetA
                    end
                end
            end #next jetA
            #    no need to do it for LH tiles, since they are implicitly done
            #    when we set NN for both jetA and jetB on the RH tiles.
        end   #next rtile
    end #next tile

    # now create the diJ (where J is i's NN) table -- remember that
    # we differ from standard normalisation here by a factor of R2
    # (corrected for at the end).
    #diJ = similar(cs.jets, _DiJ_plus_link)
    diJ = similar(cs.jets, Float64)
    NNs = similar(cs.jets, TiledJet)
    for i in eachindex(diJ)
        jetA = tiledjets[i]
        diJ[i] = _tj_diJ(jetA) # kt distance * R^2
        @inbounds NNs[i] = jetA        # our compact diJ table will not be in
        jetA.diJ_posn = i    # one-to-one corresp. with non-compact jets,
        # so set up bi-directional correspondence here.
    end

    #
    # now run the recombination loop
    history_location = length(cs.jets)
    n = length(cs.jets)
    while n > 1 #0 in C++ implementation, but the last iteration seems superfluous
        # find the minimum of the diJ on this round
        diJ_min, best = find_best(diJ, n)
        # do the recombination between A and B
        history_location += 1
        @inbounds jetA = NNs[best]
        jetB = jetA.NN

        # put the normalisation back in
        diJ_min *= invR2

        if isvalid(jetB)
            # jet-jet recombination
            # If necessary relabel A & B to ensure jetB < jetA, that way if
            # the larger of them == newtail then that ends up being jetA and
            # the new jet that is added as jetB is inserted in a position that
            # has a future!
            if jetA.id < jetB.id
                jetA, jetB = jetB, jetA;
            end

            # recombine jetA and jetB and retrieves the new index, nn
            nn = _do_ij_recombination_step!(cs, jetA.jets_index, jetB.jets_index, diJ_min)
            _tj_remove_from_tiles!(cs.tiling, jetA)
            oldB = copy(jetB)  # take a copy because we will need it...

            _tj_remove_from_tiles!(cs.tiling, jetB)
            _tj_set_jetinfo!(jetB, cs, nn, R2) # cause jetB to become _jets[nn]
            #                                  (in addition, registers the jet in the tiling)
        else
            # jet-beam recombination
            # get the hist_index
            _do_iB_recombination_step!(cs, jetA.jets_index, diJ_min)
            _tj_remove_from_tiles!(cs.tiling, jetA)
        end #isnothing(jetB)

        # first establish the set of tiles over which we are going to
        # have to run searches for updated and new nearest-neighbours --
        # basically a combination of vicinity of the tiles of the two old
        # and one new jets.
        n_near_tiles = 0
        n_near_tiles = _add_untagged_neighbours_to_tile_union(jetA.tile_index,
                                                              tile_union, n_near_tiles, tiling)
        if isvalid(jetB)
            if jetB.tile_index != jetA.tile_index
                n_near_tiles = _add_untagged_neighbours_to_tile_union(jetB.tile_index,
                                                                      tile_union, n_near_tiles, tiling)
            end
            if oldB.tile_index != jetA.tile_index && oldB.tile_index != jetB.tile_index
                n_near_tiles = _add_untagged_neighbours_to_tile_union(oldB.tile_index,
                                                                      tile_union, n_near_tiles, tiling)
            end
        end #isvalid(jetB)

        # now update our nearest neighbour info and diJ table

        # first compactify the diJ by taking the last of the diJ and copying
        # it to the position occupied by the diJ for jetA
        @inbounds NNs[n].diJ_posn = jetA.diJ_posn
        @inbounds diJ[jetA.diJ_posn] = diJ[n]
        @inbounds NNs[jetA.diJ_posn] = NNs[n]

        # then reduce size of table
        n -= 1

        # Initialise jetB's NN distance as well as updating it for
        # other particles.
        # Run over all tiles in our union
        for itile in 1:n_near_tiles
            @inbounds tile = tiling.tiles[ @inbounds tile_union[itile]] #TAKES 5μs
            @inbounds tiling.tags[tile_union[itile]] = false # reset tag, since we're done with unions

            isvalid(tile) || continue #Probably not required

            # run over all jets in the current tile
            for jetI in tile
                # see if jetI had jetA or jetB as a NN -- if so recalculate the NN
                if jetI.NN == jetA || (jetI.NN == jetB && isvalid(jetB))
                    jetI.NN_dist = R2
                    jetI.NN      = noTiledJet

                    # now go over tiles that are neighbours of I (include own tile)
                    for near_tile_index in surrounding(tile.tile_index, tiling)
                        # and then over the contents of that tile
                        for jetJ in @inbounds tiling.tiles[near_tile_index]
                            dist = _tj_dist(jetI, jetJ)
                            if dist < jetI.NN_dist && jetJ != jetI
                                jetI.NN_dist = dist
                                jetI.NN = jetJ
                            end
                        end # next jetJ
                    end # next near_tile
                    diJ[jetI.diJ_posn]  = _tj_diJ(jetI) # update diJ kt-dist
                end #jetI.NN == jetA || (jetI.NN == jetB && !isnothing(jetB))

                # check whether new jetB is closer than jetI's current NN and
                # if jetI is closer than jetB's current (evolving) nearest
                # neighbour. Where relevant update things.
                if isvalid(jetB)
                    dist = _tj_dist(jetI,jetB)
                    if dist < jetI.NN_dist
                        if jetI != jetB
                            jetI.NN_dist = dist
                            jetI.NN = jetB
                            #diJ[jetI.diJ_posn].diJ = _tj_diJ(jetI) # update diJ...
                            diJ[jetI.diJ_posn] = _tj_diJ(jetI) # update diJ...
                        end
                    end
                    if dist < jetB.NN_dist && jetI != jetB
                        jetB.NN_dist = dist
                        jetB.NN      = jetI
                    end
                end # isvalid(jetB)
            end #next jetI
        end #next itile

        # finally, register the updated kt distance for B
        if isvalid(jetB)
            @inbounds diJ[jetB.diJ_posn] = _tj_diJ(jetB)
        end
    end #next n
    inclusive_jets(cs, ptmin)
end
