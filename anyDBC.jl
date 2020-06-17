# Rewriting AnyDBC code in Julia 1.4.2
# https://github.com/VigneshN1997/AnyDBC_C_Code/blob/master/anyDBC.cpp

# Pkg.add("CSV")
# Pkg.add("Distances")
# Pkg.add("ProgressMetere")
using CSV
using Distances
using ProgressMeter

global neighbourMap = Dict{Int, Set{Int}}()
global rangeQueryPerformed = Set{Int}()
global coreList = Dict{Int, String}()
global borderList = Dict{Int, String}()
global noiseList = Set{Int}()
global coreForPointMap = Dict{Int, Set{Int}}()
global neiNoise = Set{Int}()
global clusters = Dict{Int, Set{Int}}()
global edgeNo = Dict{Int, Set{Int}}()
global edgeYes = Dict{Int, Set{Int}}()
global edgeWeak = Dict{Int, Set{Int}}()
global edgeUnknown = Dict{Int, Set{Int}}()
global edgeWeakBckup = Dict{Int, Set{Int}}()
global edgeNoBckup = Dict{Int, Set{Int}}()
global visitedNode = Dict{Int, Int}()

global usizeList = Dict{Int, Int}()
global statList = Dict{Int, Float64}()
global degList = Dict{Int, Float64}()

global popFromNoise = Vector{Int64}()
global noStatusLastIteration = Dict{Tuple{Int, Int}, Int}()

global timeElapsed = 0.0

global alpha = 512
global beta = 512
global minDist = 10
global minPts = 5


function readData(file_name::String)::Vector{Vector{Float64}}
	fp = CSV.read(file_name, delim = '\t', header=0)

	# csv format like:
	# 'id' \t 'vec'
	# vec format like:
	# 'x1, x2, x3, x4, x5, ..., xn' n=dimension

	id_list = fp[:, 1]
	vec_list = [[parse(Float64, j) for j in split(i)] for i in fp[:, 2]]
	global num_records = length(vec_list)
	global dimension = length(vec_list[1])
	println("Information about Dataset\n")
	println("Number of Records: ", num_records)
	println("Number of Dimensions: ", dimension)
	return vec_list
end

function getRandomPoints(untouchedList::Vector{Int64})::Vector{Int64}
	randomPoints = Vector{Int64}()
	if length(untouchedList) <= alpha
		append!(randomPoints, untouchedList)
	else
		ranList = rand((1:length(untouchedList)), alpha)
		append!(randomPoints, untouchedList[ranList])
	end
	return randomPoints
end

function cust_distance(x::Vector{Float64}, y::Vector{Float64})::Float64
	dist_ = Euclidean()(x, y)
	return dist_
end

function performRangeQuery(index::Int64)::Set{Int64}
	point_a = dataSet[index]
	push!(rangeQueryPerformed, index)
	neighbourPoints = Set{Int64}()
	for i in 1:num_records
		point_b = dataSet[i]
		if point_a == point_b
			continue
		end
		dist = cust_distance(point_a, point_b)
		if dist <= minDist
			push!(neighbourPoints, i)
		end
	end
	return neighbourPoints
end

function cust_append_ie(index::Int64, list::Vector{Int64})
	itr = intersect(index, list)
	if isempty(itr)
		push!(list, index)
	end
end

function cust_append_set_ie(idx1::Int64, idx2::Int64, mapList::Dict{Int64, Set{Int64}})
	itr = intersect(idx1, keys(mapList))
	if isempty(itr)
		mapList[idx1] = Set{Int64}()
	end
	push!(mapList[idx1], idx2)
end

function assignStateNei(index::Int64, listOfNeighbors::Set{Int64})::Vector{Int64}
	touchList = Vector{Int64}()
	neighbourMap[index] = listOfNeighbors
	if length(listOfNeighbors) >= minPts
		coreList[index] = "PROCESSED"
		delete!(borderList, index)

		for nei in listOfNeighbors
			cust_append_ie(nei, touchList)
			cust_append_set_ie(nei, index, coreForPointMap)

			it_rqp = intersect(nei, rangeQueryPerformed)
			if isempty(it_rqp)
				cust_append_set_ie(nei, index, neighbourMap)
				if length(neighbourMap) < minPts
					itr = intersect(nei, keys(borderList))
					if isempty(itr)
						borderList[nei] = "PROCESSED"
					end
				else
					itr_core = intersect(nei, keys(coreList))
					if isempty(itr_core)
						coreList[nei] = "UNPROCESSED"
						delete!(borderList, nei)
					end
				end
			end
		end
	else
		push!(noiseList, index)
		delete!(neiNoise, index)
		for nei in listOfNeighbors
			cust_append_ie(nei, touchList)
			itr_core = intersect(nei, keys(coreList))
			itr = intersect(nei, keys(borderList))
			if isempty(itr_core) && isempty(itr)
				push!(neiNoise, nei)
			end
			cust_append_set_ie(nei, index, neighbourMap)
		end
	end
	return touchList
end

function createPCIR(index::Int64)
	# PCIR = Primitive circle
	clusters[index] = Set{Int64}()
	push!(clusters[index], index)
end

function ddcBetPCIR(core1::Int64, core2::Int64)::Int8
	dist = cust_distance(dataSet[core1], dataSet[core2])
	if dist > 3 * minDist
		return 1
	end
	pointPCIR = neighbourMap[core1]
	push!(pointPCIR, core1)
	neiPointPCIR = neighbourMap[core2]
	push!(neiPointPCIR, core2)
	intersectionPoints = intersect(pointPCIR, neiPointPCIR)
	coreListKeys = keys(coreList)
	coreIntersectionPoints = intersect(intersectionPoints, coreListKeys)

	if !isempty(coreIntersectionPoints) || ((dist > sqrt(3) * minDist) && (length(intersectionPoints) >= minPts))
		return 0
	elseif !isempty(intersectionPoints)
		return 2
	else
		return 3
	end
end


function connComp()
	repComp = Set{Int64}()
	empty!(visitedNode)
	edgeYesKeys = keys(edgeYes)
	for keys_itr in edgeYesKeys
		visited_itr = intersect(keys_itr, keys(visitedNode))
		if isempty(visited_itr)
			DFS(keys_itr, keys_itr)
			push!(repComp, keys_itr)
		end

	end

	deleteClust = setdiff(edgeYesKeys, repComp)
	for del_itr in deleteClust
		delete!(clusters, del_itr)
	end
	empty!(edgeYes)
	count = 0
	edge_bar = Progress(length(visitedNode), dt=0.1, desc="Connectivity...")
	for visited_itr in visitedNode
		u = visited_itr[1]
		rep = visited_itr[2]
		for rep_itr in clusters[rep]
			if rep_itr == u
				continue
			end
			edgeNoItr = intersect(u, keys(edgeNo))
			if !isempty(edgeNoItr)
				delete!(edgeNo[u], rep_itr)
			end
			edgeWeakItr = intersect(u, keys(edgeWeak))
			if !isempty(edgeWeakItr)
				delete!(edgeWeak[u], rep_itr)
			end
		end
		edgeNoItr = intersect(u, keys(edgeNo))
		if !isempty(edgeNoItr)
			if isempty(edgeNo[u])
				delete!(edgeNo, u)
			end
		end
		edgeWeakItr = intersect(u, keys(edgeWeak))
		if !isempty(edgeWeakItr)
			if isempty(edgeWeak[u])
				delete!(edgeWeak, u)
			end
		end
		count += 1
		ProgressMeter.update!(edge_bar, count)
	end
end

function DFS(u::Int, rep::Int)
	visitedNode[u] = rep
	for clus_itr in clusters[u]
		push!(clusters[rep], clus_itr)
	end
	for v in edgeYes[u]
		visited_itr = intersect(v, keys(visitedNode))
		if isempty(visited_itr)
			DFS(v, rep)
		end
	end
end

function calculateStatDegree()
	w = length(clusters)
	numBorderPoints = Dict{Int64, Int64}()
	alreadyCountedPoint = Vector{Int64}()
	count = 0
	for clus_itr in clusters
		k = clus_itr[1]
		v = clus_itr[2]
		usizeList[k] = 0
		statList[k] = 0
		numBorderPoints[k] = 0
		noOfpointsInCluster = 0
		empty!(alreadyCountedPoint)
		for p in v
			for nei_itr in neighbourMap[p]
				al_itr = intersect(nei_itr, alreadyCountedPoint)
				if isempty(al_itr)
					noOfpointsInCluster += 1
					range_itr = intersect(nei_itr, rangeQueryPerformed)
					if isempty(range_itr)
						usizeList[k] += 1
					end
					border_itr = intersect(nei_itr, keys(borderList))
					if !isempty(border_itr)
						numBorderPoints[k] += 1
					end
					push!(alreadyCountedPoint, nei_itr)
				end
			end
		end
		statList[k] = (usizeList[k]/noOfpointsInCluster) + (noOfpointsInCluster/num_records)
	end
	for clus_itr in clusters
		u = clus_itr[1]
		siValue = 0
		degList[u] = 0
		edgeWeakItr = intersect(u, keys(edgeWeak))
		if !isempty(edgeWeakItr)
			for v in edgeWeak[u]
				stat_itr = intersect(v, keys(statList))
				if !isempty(stat_itr)
					degList[u] += statList[v]
					siValue += 1
				else
					delete!(edgeWeak, v)
				end
			end
			degList[u] *= w
		end
		edgeUnknownItr = intersect(u, keys(edgeUnknown))
		if !isempty(edgeUnknownItr)
			for v in edgeUnknown[u]
				stat_itr = intersect(v, keys(statList))
				if !isempty(stat_itr)
					degList[u] += statList[v]
					siValue += 1
				else
					delete!(edgeUnknown, v)
				end
			end
		end
		if numBorderPoints[u] == 0
			siValue = 0
		end
		degList[u] -= siValue
	end
end

function calculateScore()::Vector{Int64}
	scoreSet = Dict{Int64, Float64}()
	unprocessedPoints = Set{Int64}()
	range_ = Set(collect(1:num_records))

	unprocessedPoints = setdiff(range_, rangeQueryPerformed)
	unprocessedPoints1 = setdiff(unprocessedPoints, neiNoise)
	cal_bar = Progress(length(unprocessedPoints1), dt=0.05, desc="Calculating scores...")
	count = 0
	for unp_itr in unprocessedPoints1
		core_itr = intersect(unp_itr, keys(coreList))
		border_itr = intersect(unp_itr, keys(borderList))
		if !isempty(core_itr) || !isempty(border_itr)
			score = 0
			for clus_itr in clusters
				rep = clus_itr[1]
				coreListRep = clus_itr[2]
				intersectCore = intersect(coreForPointMap[unp_itr], coreListRep)
				if !isempty(intersectCore)
					score += degList[rep]
				end
			end
			score += 1/length(neighbourMap[unp_itr])
			scoreSet[unp_itr] = score
		end
		count += 1
		ProgressMeter.update!(cal_bar, count)
	end
	sorted_x = sort(collect(scoreSet), by=x->x[2])
	# sorted_x = reverse(scoreSet)
	returnList = Vector{Int64}()
	betaF = beta
	size_sorted_x = length(sorted_x)
	if size_sorted_x > 0
		if size_sorted_x < betaF
			betaF = size_sorted_x
		end
		for i in betaF:-1:1
			push!(returnList, sorted_x[i][1])
		end
	end
	return returnList
end

function dccBetPCLU(repClust1::Int64, repClust2::Int64)
	# PCLU : Primitive cluster
	noOfSubClusters = length(clusters[repClust1]) * length(clusters[repClust2])
	noCount = 0
	formedYesEdge = false
	formedWeakEdge = false
	for core1 in clusters[repClust1]
		for core2 in clusters[repClust2]
			if core1 == core2
				continue
			end
			stat = ddcBetPCIR(core1, core2)
			if stat == 0
				formedYesEdge = true
				cust_append_set_ie(repClust1, repClust2, edgeYes)
				cust_append_set_ie(repClust2, repClust1, edgeYes)
				if formedWeakEdge
					delete!(edgeWeak[repClust1], repClust2)
					if isempty(edgeWeak[repClust1])
						delete!(edgeWeak, repClust1)
					end
					delete!(edgeWeak[repClust2], repClust1)
					if isempty(edgeWeak[repClust2])
						delete!(edgeWeak, repClust2)
					end
				end
				return
			elseif stat == 1
				noCount += 1
			elseif (stat == 2) && (!formedYesEdge)
				formedWeakEdge = true
				cust_append_set_ie(repClust1, repClust2, edgeWeak)
				cust_append_set_ie(repClust2, repClust1, edgeWeak)
			end
		end
	end
	if (!formedWeakEdge && !formedYesEdge && (noCount == noOfSubClusters))
		cust_append_set_ie(repClust1, repClust2, edgeNo)
		cust_append_set_ie(repClust2, repClust1, edgeNo)
	elseif (!formedWeakEdge && !formedYesEdge && (noCount != noOfSubClusters))
		cust_append_set_ie(repClust1, repClust2, edgeUnknown)
		cust_append_set_ie(repClust2, repClust1, edgeUnknown)
	end
end

function updateStates()
	popWeakUnknown = Vector{Int64}()
	for clus_itr in clusters
		k = clus_itr[1]
		v = clus_itr[2]
		usizeList[k] = 0
		for p in v
			for x in neighbourMap[p]
				range_itr = intersect(x, rangeQueryPerformed)
				if isempty(range_itr)
					usizeList[k] += 1
				end
			end
		end
		usize_itr = intersect(k, keys(usizeList))
		if !isempty(usize_itr)
			if usizeList[k] == 0
				for edgeWeakItr in edgeWeak
					node = edgeWeakItr[1]
					neiNodes = edgeWeakItr[2]
					node_itr = intersect(k, neiNodes)
					if k == node
						for knode in neiNodes
							cust_append_set_ie(k, knode, edgeNo)
							cust_append_set_ie(knode, k, edgeNo)
						end
					elseif !isempty(node_itr)
						delete!(edgeWeak[node], k)
						if isempty(edgeWeak[node])
							push!(popWeakUnknown, node)
						end
					end
				end
				for edgeUnknownItr in edgeUnknown
					node = edgeUnknownItr[1]
					neiNodes = edgeUnknownItr[2]
					node_itr = intersect(k, neiNodes)
					if k == node
						for knode in neiNodes
							cust_append_set_ie(k, knode, edgeNo)
							cust_append_set_ie(knode, k, edgeNo)
						end
					elseif !isempty(node_itr)
						delete!(edgeUnknown[node], k)
						if isempty(edgeUnknown[node])
							push!(popWeakUnknown, node)
						end
					end
				end
				delete!(edgeWeak, k)
				delete!(edgeUnknown, k)
			end
		end
	end
	for node in popWeakUnknown
		delete!(edgeWeak, node)
		delete!(edgeUnknown, node)
	end
end

function processNoise(p::Int64)
	core_itr = intersect(p, keys(coreList))
	border_itr = intersect(p, keys(borderList))

	coreListKeys = keys(coreList)
	listOfNeighbors = neighbourMap[p]
	intersect_keys = intersect(listOfNeighbors, coreListKeys)
	if (!isempty(core_itr) || !isempty(border_itr))
		push!(popFromNoise, p)
	elseif !isempty(intersect_keys)
		borderList[p] = "PROCCESSED"
		push!(popFromNoise, p)
	else
		for nei in neighbourMap[p]
			range_itr = intersect(nei, rangeQueryPerformed)
			if !isempty(range_itr)
				listOfNeighbors = neighbourMap[nei]
			else
				listOfNeighbors = performRangeQuery(nei)
				neighbourMap[nei] = listOfNeighbors
			end
			for neiN in listOfNeighbors
				range_itr = intersect(neiN, rangeQueryPerformed)
				if isempty(range_itr)
					cust_append_set_ie(neiN, nei, neighbourMap)
					if length(neighbourMap[neiN]) >= minPts
						coreList[neiN] = "UNPROCESSED"
						delete!(borderList, neiN)
					end
				end
			end
			if length(listOfNeighbors) >= minPts
				coreList[nei] = "PROCESSED"
				borderList[p] = "PROCESSED"
				push!(popFromNoise, nei)
				push!(popFromNoise, p)
				delete!(borderList, nei)
				createPCIR(nei)
				for repCore in keys(clusters)
					if nei != repCore
						dccBetPCLU(repCore, nei)
					end
				end
				break
			end
		end
	end
end

function processOutliers()
	processNoise.(noiseList)

	for p in popFromNoise
		delete!(noiseList, p)
		delete!(neiNoise, p)
	end
	empty!(popFromNoise)
	listNei = Set{Int64}()
	for p in neiNoise
		core_itr = intersect(p, keys(coreList))
		border_itr = intersect(p, keys(borderList))
		if (!isempty(core_itr)) || (!isempty(border_itr))
			push!(popFromNoise, p)
			continue
		end
		range_itr = intersect(p, rangeQueryPerformed)
		if !isempty(range_itr)
			listNei = neighbourMap[p]
		else
			listNei = performRangeQuery(p)
			neighbourMap[p] = listNei
		end
		isCore = false
		if length(listNei) >= minPts
			isCore = true
			coreList[p] = "PROCESSED"
			delete!(borderList, p)
			push!(popFromNoise, p)
			createPCIR(p)
			for clus_itr in clusters
				repCore = clus_itr[1]
				if p != repCore
					dccBetPCLU(repCore, p)
				end
			end
		end
		for nei in listNei
			range_itr = intersect(nei, rangeQueryPerformed)
			if isempty(range_itr)
				cust_append_set_ie(nei, p, neighbourMap)
				if length(neighbourMap[nei]) >= minPts
					coreList[nei] = "UNPORCESSED"
					delete!(borderList, nei)
				end
			end
			if isCore
				borderList[nei] = "UNPORCESSED"
			end
		end
		processNoise(p)
	end
	for p in popFromNoise
		delete!(noiseList, p)
		delete!(neiNoise, p)
	end
	empty!(popFromNoise)
end


function mergeAssignNewEdge()
	empty!(edgeUnknown)
	edgeWeakBckupKeys = keys(edgeWeakBckup)
	edgeNoBckupKeys = keys(edgeNoBckup)
	eWN = union(edgeWeakBckupKeys, edgeNoBckupKeys)
	doneClusters = Set{Int64}()
	clus_bar = Progress(length(clusters), dt=0.1, desc="Merging edges...")
	count = 0
	for itr_clus1 in clusters
		k = itr_clus1[1]

		v1_1 = clusters[k]
		v1 = setdiff(v1_1, eWN)
		push!(doneClusters, k)
		for itr_clus2 in clusters
			k2 = itr_clus2[1]
			doneClusters_itr = intersect(k2, doneClusters)
			if !isempty(doneClusters_itr)
				continue
			end

			v2_1 = clusters[k2]
			v2 = setdiff(v2_1, eWN)
			weakPresent = 0
			noCount = 0
			noOfSubClusters = length(v1_1) * length(v2)
			for point in v1_1
				for neipoint in v2
					itr_point = intersect(point, keys(edgeWeak))
					itr_pointe = intersect(point, keys(edgeWeakBckup))
					if !isempty(itr_point)
						itr_neipoint = intersect(neipoint, edgeWeak[point])
						if !isempty(itr_neipoint)
							weakPresent = 1
							delete!(edgeWeak[point], itr_neipoint)
							delete!(edgeWeak[neipoint], point)
							if isempty(edgeWeak[point])
								delete!(edgeWeak, itr_point)
							end
							if isempty(edgeWeak[neipoint])
								delete!(edgeWeak, neipoint)
							end
							cust_append_set_ie(point, neipoint, edgeWeakBckup)
							cust_append_set_ie(neipoint, point, edgeWeakBckup)
						end
					elseif !isempty(itr_pointe)
						itr_neipoint = intersect(neipoint, edgeWeakBckup[point])
						if !isempty(itr_neipoint)
							weakPresent = 1
						end
					end
					itr_point = intersect(point, keys(edgeNo))
					itr_pointe = intersect(point, keys(edgeNoBckup))
					if !isempty(itr_point)
						itr_neipoint = intersect(neipoint, edgeNo[point])
						if !isempty(itr_neipoint)
							noCount += 1
							delete!(edgeNo[point], neipoint)
							delete!(edgeNo[neipoint], point)
							if isempty(edgeNo[point])
								delete!(edgeNo, point)
							end
							if isempty(edgeNo[neipoint])
								delete!(edgeNo, neipoint)
							end
							cust_append_set_ie(point, neipoint, edgeNoBckup)
							cust_append_set_ie(neipoint, point, edgeNoBckup)
						end
					elseif !isempty(itr_pointe)
						itr_neipoint = intersect(neipoint, edgeNoBckup[point])
						if !isempty(itr_neipoint)
							noCount += 1
						end
					end
				end
			end
			if weakPresent == 1
				cust_append_set_ie(k, k2, edgeWeak)
				cust_append_set_ie(k2, k, edgeWeak)
			elseif noCount == noOfSubClusters
				if (v2 == v2_1) && (v1 == v1_1)
					cust_append_set_ie(k, k2, edgeNo)
					cust_append_set_ie(k2, k, edgeNo)
					var = (k, k2)
					noStatusLastIteration[var] = noCount
				else
					noCountNew = 0
					old2 = setdiff(v2_1, v2)
					old1 = setdiff(v1_1, v1)
					for point in old2
						for neipoint in v1
							itr_point = intersect(point, keys(edgeNo))
							itr_pointe = intersect(point, keys(edgeNoBckup))
							if !isempty(itr_point)
								itr_neipoint = intersect(neipoint, edgeNo[point])
								if !isempty(itr_neipoint)
									noCountNew += 1
									delete!(edgeNo[point], neipoint)
									delete!(edgeNo[neipoint], point)
									if isempty(edgeNo[point])
										delete!(edgeNo, point)
									end
									if isempty(edgeNo[neipoint])
										delete!(edgeNo, neipoint)
									end
									cust_append_set_ie(point, neipoint, edgeNoBckup)
									cust_append_set_ie(neipoint, point, edgeNoBckup)
								end
							elseif !isempty(itr_pointe)
								itr_neipoint = intersect(neipoint, keys(edgeNoBckup))
								if !isempty(itr_neipoint)
									noCountNew += 1
								end
							end
							for p in old1
								pair1 = (point, p)
								pair2 = (p, point)
								if pair1 in keys(noStatusLastIteration)
									noCount += noStatusLastIteration[pair1]
									delete!(noStatusLastIteration, pair1)
								elseif pair2 in keys(noStatusLastIteration)
									noCount += noStatusLastIteration[pair2]
									delete!(noStatusLastIteration, pair2)
								end
							end
						end
					end
					if (noCountNew == length(v1) * length(old2))
						totalClust = noCount + noCountNew
						if totalClust == length(v1_1) * length(v2_1)
							cust_append_set_ie(k, k2, edgeNo)
							cust_append_set_ie(k2, k, edgeNo)
							var = (k, k2)
							noStatusLastIteration[var] = totalClust
						end
					end
				end
			else
				cust_append_set_ie(k, k2, edgeUnknown)
				cust_append_set_ie(k2, k, edgeUnknown)
			end
			count += 1/clusters
		end
		count += 1
		ProgressMeter.update!(clus_bar, count)
	end
	ProgressMeter.update!(clus_bar, length(clusters))
	edgeWeakKeys = keys(edgeWeak)
	edgeNoKeys = keys(edgeNo)
	edgeUnknownKeys = keys(edgeUnknown)
	clustersKeys = keys(clusters)
	removeKeys = setdiff(edgeWeakKeys, clustersKeys)
	for u in removeKeys
		delete!(edgeWeak, u)
	end
	removeKeys = setdiff(edgeNoKeys, clustersKeys)
	for u in removeKeys
		delete!(edgeNo, u)
	end
	removeKeys = setdiff(edgeUnknownKeys, clustersKeys)
	for u in removeKeys
		delete!(edgeUnknown, u)
	end
end


function stoppingCondition()::Bool
	if !isempty(edgeWeak) || !isempty(edgeUnknown)
		return true
	else
		return false
	end
end

function anyDBC()
	println("Step #1 rangeQuery on initial alpha points---")
	start_time = time()
	untouchedList = collect(1:num_records)

	chunkSize = 10
	chunkBlock = collect(1:Int64(round(num_records÷chunkSize)):num_records)[2:end-1]
	p_bar = Progress(num_records, dt=0.1, desc="Initial rangeQuery...")
	# ProgressMeter.update!(p_bar, 0.1)
	while !isempty(untouchedList)
		randomPoints = getRandomPoints(untouchedList)
		for point in randomPoints
			itr = intersect(point, untouchedList)
			if !isempty(itr)
				for j in itr
					deleteat!(untouchedList, findall(x->x == j, untouchedList))
				end
				neighbours_of_point = performRangeQuery(point)
				touch = assignStateNei(point, neighbours_of_point)
				itr_core = intersect(point, keys(coreList))
				if !isempty(itr_core)
					createPCIR(point)
				end
				for j in touch
					deleteat!(untouchedList, findall(x->x == j, untouchedList))
				end
			end
		end
		ProgressMeter.update!(p_bar, num_records - length(untouchedList))
	end
	temp_time = round(time() - start_time, digits=3)
	global timeElapsed += temp_time
	println("Step #1 rangeQuery on initial alpha points\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
	start_time = time()
	donePoints = Set()
	println("Initial count of range queries: ", length(rangeQueryPerformed), "\n")
	println("Step #2 Finding edges between clusters---")
	p_bar = Progress(length(clusters), dt=0.1, desc="Finding edges...")
	for itr_clus1 in clusters
		point = itr_clus1[1]
		push!(donePoints, point)
		for itr_clus2 in clusters
			neiPoint = itr_clus2[1]
			itr_done = intersect(neiPoint, donePoints)
			if !isempty(itr_done)
				continue
			end
			stat = ddcBetPCIR(point, neiPoint)
			if stat == 1
				cust_append_set_ie(point, neiPoint, edgeNo)
				cust_append_set_ie(neiPoint, point, edgeNo)
			elseif stat == 0
				cust_append_set_ie(point, neiPoint, edgeYes)
				cust_append_set_ie(neiPoint, point, edgeYes)
			elseif stat == 2
				cust_append_set_ie(point, neiPoint, edgeWeak)
				cust_append_set_ie(neiPoint, point, edgeWeak)
			elseif stat == 3
				cust_append_set_ie(point, neiPoint, edgeUnknown)
				cust_append_set_ie(neiPoint, point, edgeUnknown)
			end
		end
		ProgressMeter.update!(p_bar, length(donePoints))
	end
	temp_time = round(time() - start_time, digits=3)
	timeElapsed += temp_time
	println("Step #2 Finished finding edges between clusters\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
	start_time = time()
	if !isempty(edgeYes)
		println("Step #3 Finding connected components---")
		connComp()
		temp_time = round(time() - start_time, digits=3)
		timeElapsed += temp_time
		println("Step #3 Finished Finding connected components\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
		println("Step #4 Merging edges and assign new edge---")
		start_time = time()
		mergeAssignNewEdge()
		temp_time = round(time() - start_time, digits=3)
		timeElapsed += temp_time
		println("Step #4 Merge edges and assign new edge\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
		start_time = time()
	end
	iteration = 0
	while true
		println("Step #5 Checking stopping condition---")
		println("Range queries : ", length(rangeQueryPerformed), "\n")
		println("Iteration : ", iteration, "    No. of clusters : ", length(clusters), "\n")
		println("CoreList : ", length(coreList), "\n")
		println("BorderList : ", length(borderList), "\n")
		println("NoiseList : ", length(noiseList) + length(neiNoise), "\n")
		if stoppingCondition()
			temp_time = round(time() - start_time, digits=3)
			timeElapsed += temp_time
			println("Step #5 Finished checking stopping condition\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
			println("Step #6 Getting beta points for further range queries---")
			start_time = time()
			calculateStatDegree()
			betaPoints = calculateScore()
			if isempty(betaPoints)
				break
			end
			temp_time = round(time() - start_time, digits=3)
			timeElapsed += temp_time
			println("Step #6 Finished getting beta points for further range queries\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
			println("Step #7 rangeQuery on beta points, and create clusters and edges---")
			start_time = time()
			for point in betaPoints
				range_itr = intersect(point, rangeQueryPerformed)
				if !isempty(range_itr)
					listOfNeighbors = neighbourMap[point]
				else
					listOfNeighbors = performRangeQuery(point)
				end
				if length(listOfNeighbors) < minPts
					borderList[point] = "PROCESSED"
					for nei in listOfNeighbors
						range_itr = intersect(nei, rangeQueryPerformed)
						if isempty(range_itr)
							cust_append_set_ie(nei, point, neighbourMap)
						end
					end
				else
					coreList[point] = "PROCESSED"
					delete!(borderList, point)
					for nei in listOfNeighbors
						cust_append_set_ie(nei, point, coreForPointMap)
						cfpm_itr = intersect(nei, keys(coreForPointMap))
						range_itr = intersect(nei, rangeQueryPerformed)
						if isempty(range_itr)
							cust_append_set_ie(nei, point, neighbourMap)
							if length(neighbourMap[nei]) < minPts
								border_itr = intersect(nei, keys(borderList))
								if isempty(border_itr)
									range_itr = intersect(nei, rangeQueryPerformed)
									if !isempty(range_itr)
										borderList[nei] = "PROCESSED"
									else
										borderList[nei] = "UNPROCESSED"
									end
								end
							else
								core_itr = intersect(nei, keys(coreList))
								if isempty(core_itr)
									coreList[nei] = "UNPROCESSED"
									delete!(borderList, nei)
								end
							end
						end
					end
					createPCIR(point)
					for clus_itr in clusters
						repCore = clus_itr[1]
						if point != repCore
							dccBetPCLU(repCore, point)
						end
					end
				end
			end
			temp_time = round(time() - start_time, digits=3)
			timeElapsed += temp_time
			println("Step #7 rangeQuery on beta points, and create clusters and edges\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
			start_time = time()
			if !isempty(edgeYes)
				connComp()
				temp_time = round(time() - start_time, digits=3)
				timeElapsed += temp_time
				println("Step #7 Finding Connected Components\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")

				start_time = time()
				mergeAssignNewEdge()
				temp_time = round(time() - start_time, digits=3)
				timeElapsed += temp_time
				println("Step #7 Merging Connected Components\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
				start_time = time()
			end
			println("Step #8 Updating states of clusters---")
			updateStates()
			temp_time = round(time() - start_time, digits=3)
			timeElapsed += temp_time
			println("Step #8 Finished updating states of clusters if possible\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
			start_time = time()
			if !isempty(edgeYes)
				connComp()
				mergeAssignNewEdge()
				temp_time = round(time() - start_time, digits=3)
				timeElapsed += temp_time
				println("Step #8 Merge any new edges that were formed\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
				start_time = time()
			end
		else
			break
		end
		iteration += 1
	end
	processOutliers()
	temp_time = round(time() - start_time, digits=3)
	timeElapsed += temp_time
	println("Step #9 Process outliers\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
	start_time = time()
	if !isempty(edgeYes)
		connComp()
		mergeAssignNewEdge()
		temp_time = round(time() - start_time, digits=3)
		timeElapsed += temp_time
		println("Step #9 Process outliers, merge any new edges formed\n---Time used: ", temp_time, " seconds ------Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
		start_time = time()
	end
	println("Range queries : ", length(rangeQueryPerformed), "\n")
	println("CoreList : ", length(coreList), "\n")
	println("BorderList : ", length(borderList), "\n")
	println("NoiseList : ", length(noiseList) + length(neiNoise), "\n")
	range = Set()
	diff = Set()
	range = Set(collect(1:num_records))
	diff = setdiff(range, rangeQueryPerformed)
	println("Unprocessed Points : ", length(diff), "\n")
	println("No. of Clusters : ", length(clusters), "\n")
end

function main(argv)
	s_time = time()
	file_name = argv
	global dataSet = readData(file_name)
	println("\n-----AnyDBC Start------\n")
	temp_time = time() - s_time
	println("Step #0 Read data---Time used: ", round(temp_time, digits=3), " seconds ---")
	global timeElapsed += temp_time
	start_time = time()
	anyDBC()
	timeElapsed += time() - s_time
	println("---Total time elapsed: ", round(timeElapsed, digits=3), " seconds ---")
	println("------AnyDBC End-------")
	return 0
end
