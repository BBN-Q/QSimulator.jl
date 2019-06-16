"""
topleveldir
Top-level directory of QSimulator
"""
const topleveldir = dirname(@__DIR__)

"""
isbenchmarking

True if the string "bench" (case-insensitive) occurs in path to this module
"""
const isbenchmarking = occursin(r"bench"i, topleveldir)

checkbenchmarking() = isbenchmarking || error("The QSimulator module is not in a path matching r\"bench\"i.")
