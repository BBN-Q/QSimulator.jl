using Base.Test

let 
    sys = CompositeQSystem("sys")
    
    q = Qubit("q",6.)
    d = Cooling("qT1", .001, q)
    
    sys += q
    sys += d

    g = generator(sys)
end

