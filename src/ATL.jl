module ATL  # Atomic Transitions List

using Statistics, DataStructures

export get_transition_ids, get_transition, get_label, get_id, get_wavelength

abstract type Forbiddenness end
abstract type     Forbidden <: Forbiddenness  end
abstract type SemiForbidden <: Forbiddenness  end
abstract type     Permitted <: Forbiddenness  end

abstract type AbstractTransition{T <: Forbiddenness} end

struct Transition{T <: Forbiddenness} <: AbstractTransition{T}
    label::String
    lambda::Float64 # Angstrom (vacuum)
    energy_levels::NTuple{2, Float64}  # eV
end

struct Multiplet{N, T <: Forbiddenness} <: AbstractTransition{T}
    label::String
    t::NTuple{N,Transition{T}}
    function Multiplet(t::NTuple{N, Transition{T}}) where {N, T <: Forbiddenness}
        @assert N > 1
        label = t[1].label
        return new{N,T}(label, t)
    end
    function Multiplet(label::String, t::NTuple{N, Transition{T}}) where {N, T <: Forbiddenness}
        @assert N > 1
        return new{N,T}(label, t)
    end
end

struct UnidentifiedTransition <: AbstractTransition{Forbiddenness}
    label::String
    lambda::Float64 # Angstrom (vacuum)
    UnidentifiedTransition(lambda::Float64) = new("λ$(lambda)", lambda)
    UnidentifiedTransition(label::AbstractString, lambda::Float64) = new(string(label), lambda)
end

const transition_list = OrderedDict{Symbol, AbstractTransition}()

function register(t::AbstractTransition)
    global transition_list
    k = get_id(t)
    @assert !(k in keys(transition_list)) "$k is already registered"
    transition_list[k] = t
    nothing
end

# --------------------------------------------------------------------
function get_transition_ids()
    global transition_list
    return collect(keys(transition_list))
end

function get_transition(id::Symbol)
    global transition_list
    @assert id in keys(transition_list) "No transition registered with ID: $id"
    t = transition_list[id]
    @assert id == get_id(t) 
    return t
end

get_label(t::Transition) = t.label
get_label(t::Multiplet) = t.label
get_label(t::UnidentifiedTransition) = replace("l$(t.lambda)", "." => "p")
get_label(id::Symbol) = get_label(get_transition(id))

get_id(t::AbstractTransition) = Symbol(replace(get_label(t), "[" => "", "]" => "", "λ" => "_",
                                               "α" => "a", "β" => "b", "γ" => "g", "δ" => "d"))
get_id(id::Symbol) = get_id(get_transition(id))

get_wavelength(t::Transition) = t.lambda
function get_wavelength(t::Multiplet)
    id = get_id(t)
    m = mean(get_wavelength.(t.t))
    @warn "$id is a multiplet: using average wavelength: $(m)A"
    return m
end
get_wavelength(t::UnidentifiedTransition) = t.lambda
get_wavelength(id::Symbol) = get_wavelength(get_transition(id))


#=
The following data are gathered from https://www.pa.uky.edu/~peter/atomic/ using the following non-standard settings:
- Wavelength type: Vacuum
- Energy unit: eV
- Configuration (checked)
=#
#                                                                                               #  SPECTRUM  |TT|                       CONFIGURATION                       |        TERM         |   JJ    |    LEVEL_ENERGY_eV       |
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Transition{Permitted}(    "Lyδ"          ,   949.74298 , (  0.000000,  13.054507)))    #  H I       |E1|                           1*-5*                           |         1-5         | 1/2-*   |    0.000000 -   13.054507|
register(Transition{Permitted}(    "Lyγ"          ,   972.53674 , (  0.000000,  12.748544)))    #  H I       |E1|                           1*-4*                           |         1-4         | 1/2-*   |    0.000000 -   12.748544|
register(Transition{Permitted}(    "Lyβ"          ,  1025.7222  , (  0.000000,  12.087510)))    #  H I       |E1|                           1*-3*                           |         1-3         | 1/2-*   |    0.000000 -   12.087510|
register(Transition{Permitted}(    "Lyα"          ,  1215.6700  , (  0.000000,  10.198834)))    #  H I       |E1|                           1*-2*                           |         1-2         | 1/2-*   |    0.000000 -   10.198834|
register(Transition{Permitted}(    "Hd"           ,  4102.892   , ( 10.198834,  13.220709)))    #  H I       |E1|                           2*-6*                           |         2-6         |   *-*   |   10.198834 -   13.220709|
register(Transition{Permitted}(    "Hγ"           ,  4341.684   , ( 10.198834,  13.054507)))    #  H I       |E1|                           2*-5*                           |         2-5         |   *-*   |   10.198834 -   13.054507|
register(Transition{Permitted}(    "Hβ"           ,  4862.683   , ( 10.198834,  12.748544)))    #  H I       |E1|                           2*-4*                           |         2-4         |   *-*   |   10.198834 -   12.748544|
register(Transition{Permitted}(    "Hα"           ,  6564.61    , ( 10.198834,  12.087510)))    #  H I       |E1|                           2*-3*                           |         2-3         |   *-*   |   10.198834 -   12.087510|
register(Transition{Permitted}(    "Paβ"          , 12821.59    , ( 12.087510,  13.054507)))    #  H I       |E1|                           3*-5*                           |         3-5         |   *-*   |   12.087510 -   13.054507|
register(Transition{Permitted}(    "Paα"          , 18756.13    , ( 12.087510,  12.748544)))    #  H I       |E1|                           3*-4*                           |         3-4         |   *-*   |   12.087510 -   12.748544|
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Transition{Permitted}(    "HeIIλ1215"    ,  1215.134   , ( 40.813419,  51.016783)))    #  He II     |E1|                           2*-4*                           |         2-4         |   *-*   |   40.813419 -   51.016783|
register(Transition{Permitted}(    "HeIIλ1640"    ,  1640.42    , ( 40.813419,  48.371511)))    #  He II     |E1|                           2*-3*                           |         2-3         |   *-*   |   40.813419 -   48.371511|
register(Multiplet((
         Transition{Permitted}(    "HeIλ4471"     ,  4472.72514 , ( 20.964108,  23.736115)),    #  He I      |E1|                        1s.2p-1s.4d                        |       3Po-3D        |   2-1   |   20.964108 -   23.736115|
         Transition{Permitted}(    "HeIλ4471"     ,  4472.72884 , ( 20.964108,  23.736112)),    #  He I      |E1|                        1s.2p-1s.4d                        |       3Po-3D        |   2-2   |   20.964108 -   23.736112|
         Transition{Permitted}(    "HeIλ4471"     ,  4472.72909 , ( 20.964108,  23.736112)),    #  He I      |E1|                        1s.2p-1s.4d                        |       3Po-3D        |   2-3   |   20.964108 -   23.736112|
         Transition{Permitted}(    "HeIλ4471"     ,  4472.74043 , ( 20.964117,  23.736115)),    #  He I      |E1|                        1s.2p-1s.4d                        |       3Po-3D        |   1-1   |   20.964117 -   23.736115|
         Transition{Permitted}(    "HeIλ4471"     ,  4472.74413 , ( 20.964117,  23.736112)),    #  He I      |E1|                        1s.2p-1s.4d                        |       3Po-3D        |   1-2   |   20.964117 -   23.736112|
         Transition{Permitted}(    "HeIλ4471"     ,  4472.93808 , ( 20.964240,  23.736115)))))  #  He I      |E1|                        1s.2p-1s.4d                        |       3Po-3D        |   0-1   |   20.964240 -   23.736115|
register(Transition{Permitted}(    "HeIIλ4686"    ,  4687.02    , ( 48.371511,  51.016783)))    #  He II     |E1|                           3*-4*                           |         3-4         |   *-*   |   48.371511 -   51.016783|
register(Multiplet((
         Transition{Permitted}(    "HeIλ5876"     ,  5877.227156, ( 20.964108,  23.073678)),    #  He I      |E1|                        1s.2p-1s.3d                        |       3Po-3D        |   2-1   |   20.964108 -   23.073678|
         Transition{Permitted}(    "HeIλ5876"     ,  5877.242417, ( 20.964108,  23.073673)),    #  He I      |E1|                        1s.2p-1s.3d                        |       3Po-3D        |   2-2   |   20.964108 -   23.073673|
         Transition{Permitted}(    "HeIλ5876"     ,  5877.243294, ( 20.964108,  23.073673)),    #  He I      |E1|                        1s.2p-1s.3d                        |       3Po-3D        |   2-3   |   20.964108 -   23.073673|
         Transition{Permitted}(    "HeIλ5876"     ,  5877.253555, ( 20.964117,  23.073678)),    #  He I      |E1|                        1s.2p-1s.3d                        |       3Po-3D        |   1-1   |   20.964117 -   23.073678|
         Transition{Permitted}(    "HeIλ5876"     ,  5877.268816, ( 20.964117,  23.073673)),    #  He I      |E1|                        1s.2p-1s.3d                        |       3Po-3D        |   1-2   |   20.964117 -   23.073673|
         Transition{Permitted}(    "HeIλ5876"     ,  5877.594821, ( 20.964240,  23.073678)))))  #  He I      |E1|                        1s.2p-1s.3d                        |       3Po-3D        |   0-1   |   20.964240 -   23.073678|
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Multiplet((
         Transition{Permitted}(    "CIIλ1335"     ,  1334.5323  , (  0.000000,   9.290464)),    #  C II      |E1|                       2s2.2p-2s.2p2                       |       2Po-2D        | 1/2-3/2 |    0.000000 -    9.290464|
         Transition{Permitted}(    "CIIλ1335"     ,  1335.6627  , (  0.007863,   9.290464)),    #  C II      |E1|                       2s2.2p-2s.2p2                       |       2Po-2D        | 3/2-3/2 |    0.007863 -    9.290464|
         Transition{Permitted}(    "CIIλ1335"     ,  1335.7077  , (  0.007863,   9.290152)))))  #  C II      |E1|                       2s2.2p-2s.2p2                       |       2Po-2D        | 3/2-5/2 |    0.007863 -    9.290152|
register(Multiplet((
         Transition{SemiForbidden}("CII]λ2326"    ,  2324.21    , (  0.000000,   5.334459)),    #  C II]     |E1|                       2s2.2p-2s.2p2                       |       2Po-4P        | 1/2-3/2 |    0.000000 -    5.334459|
         Transition{SemiForbidden}("CII]λ2326"    ,  2325.40    , (  0.000000,   5.331732)),    #  C II]     |E1|                       2s2.2p-2s.2p2                       |       2Po-4P        | 1/2-1/2 |    0.000000 -    5.331732|
         Transition{SemiForbidden}("CII]λ2326"    ,  2326.11    , (  0.007863,   5.337968)),    #  C II]     |E1|                       2s2.2p-2s.2p2                       |       2Po-4P        | 3/2-5/2 |    0.007863 -    5.337968|
         Transition{SemiForbidden}("CII]λ2326"    ,  2327.64    , (  0.007863,   5.334459)),    #  C II]     |E1|                       2s2.2p-2s.2p2                       |       2Po-4P        | 3/2-3/2 |    0.007863 -    5.334459|
         Transition{SemiForbidden}("CII]λ2326"    ,  2328.84    , (  0.007863,   5.331732)))))  #  C II]     |E1|                       2s2.2p-2s.2p2                       |       2Po-4P        | 3/2-1/2 |    0.007863 -    5.331732|
register(Multiplet((
         Transition{Permitted}(    "CIVλ1549"     ,  1548.203   , (  0.000000,   8.008266)),    #  C IV      |E1|                       1s2.2s-1s2.2p                       |        2S-2Po       | 1/2-3/2 |    0.000000 -    8.008266|
         Transition{Permitted}(    "CIVλ1549"     ,  1550.777   , (  0.000000,   7.994975)))))  #  C IV      |E1|                       1s2.2s-1s2.2p                       |        2S-2Po       | 1/2-1/2 |    0.000000 -    7.994975|
register(Transition{SemiForbidden}("CIIIλ1909"    ,  1908.53    , ( 44.275974,  50.772292)))    #  C III     |E1|                        2s.6p-2p.5p                        |       1Po-1D        |   1-2   |   44.275974 -   50.772292|  Transitions from auto-ionizing levels. TODO: is it correct?
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Multiplet((
         Transition{Permitted}(    "AlIIIλ1858"   ,   1854.716  , (  0.000000,   6.684809)),    #  Al III    |E1|                           3s-3p                           |        2S-2Po       | 1/2-3/2 |    0.000000 -    6.684809|
         Transition{Permitted}(    "AlIIIλ1858"   ,   1862.790  , (  0.000000,   6.655838)))))  #  Al III    |E1|                           3s-3p                           |        2S-2Po       | 1/2-1/2 |    0.000000 -    6.655838|
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Transition{Forbidden}(    "[NII]λ6549"   ,  6549.85    , (  0.006034,   1.898966)))    # [N II]     |M1|                      2s2.2p2-2s2.2p2                      |        3P-1D        |   1-2   |    0.006034 -    1.898966|
register(Transition{Forbidden}(    "[NII]λ6583"   ,  6585.28    , (  0.016217,   1.898966)))    # [N II]     |M1|                      2s2.2p2-2s2.2p2                      |        3P-1D        |   2-2   |    0.016217 -    1.898966|
register(Multiplet((
         Transition{Permitted}(    "NVλ1241"      ,  1238.821   , (  0.000000,  10.008244)),    #  N V       |E1|                       1s2.2s-1s2.2p                       |        2S-2Po       | 1/2-3/2 |    0.000000 -   10.008244|
         Transition{Permitted}(    "NVλ1241"      ,  1242.804   , (  0.000000,   9.976169)))))  #  N V       |E1|                       1s2.2s-1s2.2p                       |        2S-2Po       | 1/2-1/2 |    0.000000 -    9.976169|
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Transition{Forbidden}(    "[NeV]λ3345"   ,  3346.4     , (  0.050986,   3.755979)))    # [Ne V]     |M1|                      2s2.2p2-2s2.2p2                      |        3P-1D        |   1-2   |    0.050986 -    3.755979|
register(Transition{Forbidden}(    "[NeV]λ3426"   ,  3426.5     , (  0.137557,   3.755979)))    # [Ne V]     |M1|                      2s2.2p2-2s2.2p2                      |        3P-1D        |   2-2   |    0.137557 -    3.755979|
register(Transition{Forbidden}(    "[NeIII]λ3869" ,  3870.16    , (  0.000000,   3.203592)))    # [Ne III]   |M1|                      2s2.2p4-2s2.2p4                      |        3P-1D        |   2-2   |    0.000000 -    3.203592|
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Transition{Forbidden}(    "[OI]λ6300"    ,  6302.046   , (  0.000000,   1.967365)))    # [O I]      |M1|                      2s2.2p4-2s2.2p4                      |        3P-1D        |   2-2   |    0.000000 -    1.967365|
register(Transition{Forbidden}(    "[OI]λ6364"    ,  6365.536   , (  0.019622,   1.967365)))    # [O I]      |M1|                      2s2.2p4-2s2.2p4                      |        3P-1D        |   1-2   |    0.019622 -    1.967365|
register(Multiplet((
         Transition{Forbidden}(    "[OII]λ3727"   ,  3727.092   , (  0.000000,   3.326568)),    # [O II]     |E2|                      2s2.2p3-2s2.2p3                      |       4So-2Do       | 3/2-3/2 |    0.000000 -    3.326568|
         Transition{Forbidden}(    "[OII]λ3727"   ,  3729.875   , (  0.000000,   3.324086)))))  # [O II]     |E2|                      2s2.2p3-2s2.2p3                      |       4So-2Do       | 3/2-5/2 |    0.000000 -    3.324086|
register(Multiplet((
         Transition{Permitted}(    "OIλ1306"      ,  1302.16848 , (  0.000000,   9.521367)),    # O I        |E1|                   2s2.2p4-2s2.2p3.(4So).3s                |        3P-3So       |   2-1   |    0.000000 -    9.521367|
         Transition{Permitted}(    "OIλ1306"      ,  1304.85763 , (  0.019622,   9.521367)),    # O I        |E1|                   2s2.2p4-2s2.2p3.(4So).3s                |        3P-3So       |   1-1   |    0.019622 -    9.521367|
         Transition{Permitted}(    "OIλ1306"      ,  1306.02861 , (  0.028142,   9.521367)))))  # O I        |E1|                   2s2.2p4-2s2.2p3.(4So).3s                |        3P-3So       |   0-1   |    0.028142 -    9.521367|
register(Multiplet((
         Transition{SemiForbidden}("OIII]λ1664"   ,  1660.8092  , (  0.014032,   7.479324)),    # O III]     |E1|                       2s2.2p2-2s.2p3                      |        3P-5So       |   1-2   |    0.014032 -    7.479324|
         Transition{SemiForbidden}("OIII]λ1664"   ,  1666.1497  , (  0.037961,   7.479324)))))  # O III]     |E1|                       2s2.2p2-2s.2p3                      |        3P-5So       |   2-2   |    0.037961 -    7.479324|
register(Transition{Forbidden}(    "[OIII]λ4363"  ,  4364.436   , (  2.513566,   5.354351)))    # [O III]    |E2|                      2s2.2p2-2s2.2p2                      |        1D-1S        |   2-0   |    2.513566 -    5.354351|
register(Transition{Forbidden}(    "[OIII]λ4932"  ,  4932.603   , (  0.000000,   2.513566)))    # [O III]    |E2|                      2s2.2p2-2s2.2p2                      |        3P-1D        |   0-2   |    0.000000 -    2.513566|
register(Transition{Forbidden}(    "[OIII]λ4959"  ,  4960.295   , (  0.014032,   2.513566)))    # [O III]    |M1|                      2s2.2p2-2s2.2p2                      |        3P-1D        |   1-2   |    0.014032 -    2.513566|
register(Transition{Forbidden}(    "[OIII]λ5007"  ,  5008.240   , (  0.037961,   2.513566)))    # [O III]    |M1|                      2s2.2p2-2s2.2p2                      |        3P-1D        |   2-2   |    0.037961 -    2.513566|
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Multiplet((
         Transition{Permitted}(    "MgIIλ2798"    ,  2796.352   , (  0.000000,   4.433786)),    #  Mg II     |E1|                           3s-3p                           |        2S-2Po       | 1/2-3/2 |    0.000000 -    4.433786|
         Transition{Permitted}(    "MgIIλ2798"    ,  2803.531   , (  0.000000,   4.422432)))))  #  Mg II     |E1|                           3s-3p                           |        2S-2Po       | 1/2-1/2 |    0.000000 -    4.422432|
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Multiplet((
         Transition{Permitted}(    "SiIVλ1400"    ,  1393.7546  , (  0.000000,   8.895701)),    # Si IV      |E1|                           3s-3p                           |        2S-2Po       | 1/2-3/2 |    0.000000 -    8.895701|
         Transition{Permitted}(    "SiIVλ1400"    ,  1402.7697  , (  0.000000,   8.838532)))))  # Si IV      |E1|                           3s-3p                           |        2S-2Po       | 1/2-1/2 |    0.000000 -    8.838532|
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
register(Transition{Forbidden}(    "SIIλ6716"     ,  6718.29    , (  0.000000,   1.845472)))    # [S II]     |E2|                      3s2.3p3-3s2.3p3                      |       4So-2Do       | 3/2-5/2 |    0.000000 -    1.845472|
register(Transition{Forbidden}(    "SIIλ6731"     ,  6732.67    , (  0.000000,   1.841531)))    # [S II]     |E2|                      3s2.3p3-3s2.3p3                      |       4So-2Do       | 3/2-3/2 |    0.000000 -    1.841531|

end
