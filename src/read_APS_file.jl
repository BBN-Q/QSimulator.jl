using HDF5

function read_APS_file(filename, startIdx=0)

	#APS bit masks
	START_MINILL_MASK = 2^15
	END_MINILL_MASK = 2^14
	TA_PAIR_MASK = 2^12
	REPEAT_MASK = 2^10 - 1
	MAX_WAVEFORM_VALUE = 2^13 - 1
	ADDRESS_UNIT = 4

	chanStrs = ["chan_1", "chan_2", "chan_3", "chan_4"]
	mrkStrs = ["ch1m1", "ch2m1", "ch3m1", "ch4m1"]
	AWGData = Dict(chan => Vector{Float64}[] for chan in chanStrs)

	h5open(filename, "r") do f
		for (ct, chanStr) in enumerate(chanStrs)
			if read(attrs(f[chanStr])["isIQMode"]) != 0
				curLLData = f[chanStrs[2*fld(ct-1, 2) + 1]]["linkListData"]
			else
				curLLData = f[chanStr]["linkListData"]
			end

			# pull out LL data
			addr = read(curLLData["addr"]) * ADDRESS_UNIT
			count = (1 .+ read(curLLData["count"])) * ADDRESS_UNIT
			repeat = read(curLLData["repeat"])
			trigger1 = read(curLLData["trigger1"]) * ADDRESS_UNIT
			trigger2 = read(curLLData["trigger2"]) * ADDRESS_UNIT
			numEntries = read(attrs(curLLData)["length"])[1]

			# pull out and scale waveform data
			wfLib = (1./MAX_WAVEFORM_VALUE)*read(f[chanStr]["waveformLib"])

			# initialize storage
			AWGData[mrkStrs[ct]] = Vector{Float64}[]

			# loop over LL entries
			for entryct = 1:numEntries
				# at the start of a miniLL, create an empty array
				if repeat[entryct] & START_MINILL_MASK != 0
					curWF = Float64[]
					triggerDelays = Int[]
				end

				# record the trigger delays
				if ct % 2 == 1
					if trigger1[entryct] > 0
						push!(triggerDelays, length(curWF) + trigger1[entryct])
					end
				else
					if trigger2[entryct] > 0
						push!(triggerDelays, length(curWF) + trigger2[entryct])
					end
				end

				# if it is a TA pair or a regular pulse
				curRepeat = (repeat[entryct] & REPEAT_MASK) + 1
				if repeat[entryct] & TA_PAIR_MASK != 0
					curWF = [curWF, wfLib[1+addr[entryct]]*ones(count[entryct] * curRepeat)]
				else
					try
						curWF = [curWF, repmat(wfLib[1+addr[entryct]:addr[entryct]+count[entryct]], curRepeat, 1)]
					catch
						print("addr: ")
						println(int(addr[entryct]))
						print("count: ")
						println(int(count[entryct]))
					end
				end

				# add the trigger pulses and push on the curWF
				if repeat[entryct] & END_MINILL_MASK != 0
					# println("Size of curWF: $(size(curWF))")
					push!(AWGData[chanStr], vec(curWF))
					triggerSeq = zeros(Uint8, length(curWF))
					triggerSeq[triggerDelays] = 1
					push!(AWGData[mrkStrs[ct]], triggerSeq)
				end
			end
		end
	end

	#Strip-off the starting zero points
	if startIdx > 0
		for seqs in values(AWGData)
			for seq in seqs
				deleteat!(seq, 1:startIdx)
			end
		end
	end

	return AWGData
end
