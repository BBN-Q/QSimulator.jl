using QSimulator

using Base.Test

r = Resonator("silly",  5.5, 3)
@test QSimulator.raising(r) * QSimulator.lowering(r) == QSimulator.number(r)
