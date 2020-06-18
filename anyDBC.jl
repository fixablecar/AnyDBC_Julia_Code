# Rewriting AnyDBC code in Julia 1.4.2
# https://github.com/VigneshN1997/AnyDBC_C_Code/blob/master/anyDBC.cpp

using CSV
using Distances
using ProgressMeter
using ThreadPools

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
	for i in setdiff(1:num_records, index)
		point_b = dataSet[i]
		dist = cust_distance(point_a, point_b)
		if dist <= minDist
			push!(neighbourPoints, i)
		end
	end
	return neighbourPoints
end

function pushList(index::Int64, list::Vector{Int64})
	if index ∉ list
		push!(list, index)
	end
end

function pushDict(idx1::Int64, idx2::Int64, mapList::Dict{Int64, Set{Int64}})
	if !haskey(mapList, idx1)
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
			pushList(nei, touchList)
			pushDict(nei, index, coreForPointMap)
			if nei ∉ rangeQueryPerformed
				pushDict(nei, index, neighbourMap)
				if length(neighbourMap) < minPts
					if !haskey(borderList, nei)
						borderList[nei] = "PROCESSED"
					end
				else
					if !haskey(coreList, nei)
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
			pushList(nei, touchList)
			if !haskey(coreList, nei) && !haskey(borderList, nei)
				push!(neiNoise, nei)
			end
			pushDict(nei, index, neighbourMap)
		end
	end
	return touchList
end

function createPCIR(index::Int64)
	# PCIR = Primitive circle
	pushDict(index, index, clusters)
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
	edge_bar = Progress(length(edgeYesKeys), dt=0.1, desc="Connectivity...")
	count = 0
	for keys_itr in edgeYesKeys
		if !haskey(visitedNode, keys_itr)
			DFS(keys_itr, keys_itr)
			push!(repComp, keys_itr)
		end
		count += 1
		ProgressMeter.update!(edge_bar, count)
	end
	deleteClust = setdiff(edgeYesKeys, repComp)
	for del_itr in deleteClust
		delete!(clusters, del_itr)
	end
	empty!(edgeYes)
	for visited_itr in visitedNode
		u = visited_itr[1]
		rep = visited_itr[2]
		rep_range = setdiff(clusters[rep], u)
		for rep_itr in rep_range
			if haskey(edgeNo, u)
				delete!(edgeNo[u], rep_itr)
			end
			if haskey(edgeWeak, u)
				delete!(edgeWeak[u], rep_itr)
			end
		end
		if haskey(edgeNo, u)
			if isempty(edgeNo[u])
				delete!(edgeNo, u)
			end
		end
		if haskey(edgeWeak, u)
			if isempty(edgeWeak[u])
				delete!(edgeWeak, u)
			end
		end
	end
end

function DFS(u::Int, rep::Int)
	visitedNode[u] = rep
	for v in clusters[u]
		push!(clusters[rep], v)
	end
	for v in edgeYes[u]
		if !haskey(visitedNode, v)
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
			for x in neighbourMap[p]
				if x ∉ alreadyCountedPoint
					noOfpointsInCluster += 1
					if x ∉ rangeQueryPerformed
						usizeList[k] += 1
					end
					if haskey(borderList, x)
						numBorderPoints[k] += 1
					end
					append!(alreadyCountedPoint, x)
				end
			end
		end
		statList[k] = (usizeList[k]/noOfpointsInCluster) + (noOfpointsInCluster/num_records)
	end
	for clus_itr in clusters
		u = clus_itr[1]
		siValue = 0
		degList[u] = 0
		if haskey(edgeWeak, u)
			for v in edgeWeak[u]
				if haskey(statList, v)
					degList[u] += statList[v]
					siValue += 1
				else
					delete!(edgeWeak, v)
				end
			end
			degList[u] *= w
		end
		if haskey(edgeUnknown, u)
			for v in edgeUnknown[u]
				if haskey(statList, v)
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
	unprocessedPoints = setdiff(1:num_records, rangeQueryPerformed)
	unprocessedPoints1 = setdiff(unprocessedPoints, neiNoise)
	cal_bar = Progress(length(unprocessedPoints1), dt=0.05, desc="Calculating scores...")
	count = 0
	for unp_itr in unprocessedPoints1
		if haskey(coreList, unp_itr) || haskey(borderList, unp_itr)
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
	sorted_x = sort(unique(scoreSet), by=x->x[2])
	returnList = Vector{Int64}()
	betaF = beta
	size_sorted_x = length(sorted_x)
	if size_sorted_x > 0
		if size_sorted_x < betaF
			betaF = size_sorted_x
		end
		for i in betaF:-1:1
			append!(returnList, sorted_x[i][1])
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
		core2_range = setdiff(clusters[repClust2], core1)
		for core2 in core2_range
			stat = ddcBetPCIR(core1, core2)
			if stat == 0
				formedYesEdge = true
				pushDict(repClust1, repClust2, edgeYes)
				pushDict(repClust2, repClust1, edgeYes)
				if formedWeakEdge
					clearPointPairs((repClust1, repClust2), edgeWeak)
				end
				return
			elseif stat == 1
				noCount += 1
			elseif (stat == 2) && (!formedYesEdge)
				formedWeakEdge = true
				pushDict(repClust1, repClust2, edgeWeak)
				pushDict(repClust2, repClust1, edgeWeak)
			end
		end
	end
	if (!formedWeakEdge && !formedYesEdge && (noCount == noOfSubClusters))
		pushDict(repClust1, repClust2, edgeNo)
		pushDict(repClust2, repClust1, edgeNo)
	elseif (!formedWeakEdge && !formedYesEdge && (noCount != noOfSubClusters))
		pushDict(repClust1, repClust2, edgeUnknown)
		pushDict(repClust2, repClust1, edgeUnknown)
	end
end

function updateStates()
	popWeakUnknown = Vector{Int64}()
	for k in keys(clusters)
		v = clusters[k]
		usizeList[k] = 0
		for p in v
			for x in neighbourMap[p]
				if x ∉ rangeQueryPerformed
					usizeList[k] += 1
				end
			end
		end
		if haskey(usizeList, k)
			if usizeList[k] == 0
				for node in keys(edgeWeak)
					neiNodes = edgeWeak[node]
					if k == node
						for knode in neiNodes
							pushDict(k, knode, edgeNo)
							pushDict(knode, k, edgeNo)
						end
					elseif k ∈ neiNodes
						delete!(edgeWeak[node], k)
						if isempty(edgeWeak[node])
							append!(popWeakUnknown, node)
						end
					end
				end
				for node in keys(edgeUnknown)
					neiNodes = edgeUnknown[node]
					if k == node
						for knode in neiNodes
							pushDict(k, knode, edgeNo)
							pushDict(knode, k, edgeNo)
						end
					elseif k ∈ neiNodes
						delete!(edgeUnknown[node], k)
						if isempty(edgeUnknown[node])
							append!(popWeakUnknown, node)
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
	coreListKeys = keys(coreList)
	listOfNeighbors = neighbourMap[p]
	intersect_keys = intersect(listOfNeighbors, coreListKeys)
	if haskey(coreList, p) || haskey(borderList, p)
		append!(popFromNoise, p)
	elseif !isempty(intersect_keys)
		borderList[p] = "PROCCESSED"
		append!(popFromNoise, p)
	else
		for nei in neighbourMap[p]
			if nei ∈ rangeQueryPerformed
				listOfNeighbors = neighbourMap[nei]
			else
				listOfNeighbors = performRangeQuery(nei)
				neighbourMap[nei] = listOfNeighbors
			end
			for neiN in listOfNeighbors
				if neiN ∉ rangeQueryPerformed
					pushDict(neiN, nei, neighbourMap)
					if length(neighbourMap[neiN]) >= minPts
						coreList[neiN] = "UNPROCESSED"
						delete!(borderList, neiN)
					end
				end
			end
			if length(listOfNeighbors) >= minPts
				coreList[nei] = "PROCESSED"
				borderList[p] = "PROCESSED"
				append!(popFromNoise, nei)
				append!(popFromNoise, p)
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
		if haskey(coreList, p) || haskey(borderList, p)
			append!(popFromNoise, p)
		else
			if p ∈ rangeQueryPerformed
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
				append!(popFromNoise, p)
				createPCIR(p)
				for repCore in keys(clusters)
					if p != repCore
						dccBetPCLU(repCore, p)
					end
				end
			end
			for nei in listNei
				if nei ∉ rangeQueryPerformed
					pushDict(nei, p, neighbourMap)
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
	end
	for p in popFromNoise
		delete!(noiseList, p)
		delete!(neiNoise, p)
	end
	empty!(popFromNoise)
end

function clearPointPairs(pt_pair, edge)
	point, neipoint = pt_pair
	delete!(edge[point], neipoint)
	if isempty(edge[point])
		delete!(edge, point)
	end
	delete!(edge[neipoint], point)
	if isempty(edge[neipoint])
		delete!(edge, neipoint)
	end
end

function mergeAssignNewEdge()
	itr_del = Vector{Tuple{Int64, Int64}}()
	empty!(edgeUnknown)
	eWN = union(keys(edgeWeakBckup), keys(edgeNoBckup))
	doneClusters = Set{Int64}()
	clustersKeys = keys(clusters)
	for k in clustersKeys
		clus_bar = Progress(length(clusters), desc="Merging edges...#$k ")
		count = 0
		ProgressMeter.update!(clus_bar, count)
		v1_1 = clusters[k]
		v1 = setdiff(v1_1, eWN)
		push!(doneClusters, k)
		for k2 in clustersKeys
			if k2 ∉ doneClusters
				v2_1 = clusters[k2]
				v2 = setdiff(v2_1, eWN)
				weakPresent = 0
				noCount = 0
				noOfSubClusters = length(v1_1) * length(v2)
				for point in v1_1
					for neipoint in v2
						if haskey(edgeWeak, point)
							if neipoint ∈ edgeWeak[point]
								weakPresent = 1
								clearPointPairs((point, neipoint), edgeWeak)
								pushDict(point, neipoint, edgeWeakBckup)
								pushDict(neipoint, point, edgeWeakBckup)
							end
						elseif haskey(edgeWeakBckup, point)
							if neipoint ∈ edgeWeakBckup[point]
								weakPresent = 1
							end
						end
						if haskey(edgeNo, point)
							if neipoint ∈ edgeNo[point]
								noCount += 1
								clearPointPairs((point, neipoint), edgeNo)
								pushDict(point, neipoint, edgeNoBckup)
								pushDict(neipoint, point, edgeNoBckup)
							end
						elseif haskey(edgeNoBckup, point)
							if neipoint ∈ edgeNoBckup[point]
								noCount += 1
							end
						end
					end
				end
				if weakPresent == 1
					pushDict(k, k2, edgeWeak)
					pushDict(k2, k, edgeWeak)
				elseif noCount == noOfSubClusters
					if (v2 == v2_1) && (v1 == v1_1)
						pushDict(k, k2, edgeNo)
						pushDict(k2, k, edgeNo)
						var = (k, k2)
						noStatusLastIteration[var] = noCount
					else
						noCountNew = 0
						old2 = setdiff(v2_1, v2)
						old1 = setdiff(v1_1, v1)
						for point in old2
							for neipoint in v1
								if haskey(edgeNo, point)
									if neipoint ∈ edgeNo[point]
										noCountNew += 1
										clearPointPairs((point, neipoint), edgeNo)
										pushDict(point, neipoint, edgeNoBckup)
										pushDict(neipoint, point, edgeNoBckup)
									end
								elseif haskey(edgeNoBckup, point)
									if haskey(edgeNoBckup, neipoint)
										noCountNew += 1
									end
								end
								for p in old1
									pair1 = (point, p)
									pair2 = (p, point)
									if pair1 in keys(noStatusLastIteration)
										noCount += noStatusLastIteration[pair1]
										push!(itr_del, pair1)
									elseif pair2 in keys(noStatusLastIteration)
										noCount += noStatusLastIteration[pair2]
										push!(itr_del, pair2)
									end
								end
							end
						end
						for ptx in itr_del
							delete!(noStatusLastIteration, ptx)
						end
						if noCountNew == (length(v1) * length(old2))
							totalClust = noCount + noCountNew
							if totalClust == (length(v1_1) * length(v2_1))
								pushDict(k, k2, edgeNo)
								pushDict(k2, k, edgeNo)
								var = (k, k2)
								noStatusLastIteration[var] = totalClust
							end
						end
					end
				else
					pushDict(k, k2, edgeUnknown)
					pushDict(k2, k, edgeUnknown)
				end
			end
			count += 1
			ProgressMeter.update!(clus_bar, count)
		end
	end
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
			if point ∈ untouchedList
				setdiff!(untouchedList, point)
				neighbours_of_point = performRangeQuery(point)
				touch = assignStateNei(point, neighbours_of_point)
				if haskey(coreList, point)
					createPCIR(point)
				end
				setdiff!(untouchedList, touch)
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
	clusters_keys = keys(clusters)
	for point in clusters_keys
		push!(donePoints, point)
		for neiPoint in clusters_keys
			if neiPoint ∉ donePoints
				stat = ddcBetPCIR(point, neiPoint)
				if stat == 1
					pushDict(point, neiPoint, edgeNo)
					pushDict(neiPoint, point, edgeNo)
				elseif stat == 0
					pushDict(point, neiPoint, edgeYes)
					pushDict(neiPoint, point, edgeYes)
				elseif stat == 2
					pushDict(point, neiPoint, edgeWeak)
					pushDict(neiPoint, point, edgeWeak)
				elseif stat == 3
					pushDict(point, neiPoint, edgeUnknown)
					pushDict(neiPoint, point, edgeUnknown)
				end
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
				if point ∈ rangeQueryPerformed
					listOfNeighbors = neighbourMap[point]
				else
					listOfNeighbors = performRangeQuery(point)
				end
				if length(listOfNeighbors) < minPts
					borderList[point] = "PROCESSED"
					for nei in listOfNeighbors
						if nei ∉ rangeQueryPerformed
							pushDict(nei, point, neighbourMap)
						end
					end
				else
					coreList[point] = "PROCESSED"
					delete!(borderList, point)
					for nei in listOfNeighbors
						pushDict(nei, point, coreForPointMap)
						if nei ∉ rangeQueryPerformed
							pushDict(nei, point, neighbourMap)
							if length(neighbourMap[nei]) < minPts
								if !haskey(borderList, nei)
									if nei ∈ rangeQueryPerformed
										borderList[nei] = "PROCESSED"
									else
										borderList[nei] = "UNPROCESSED"
									end
								end
							else
								if !haskey(coreList, nei)
									coreList[nei] = "UNPROCESSED"
									delete!(borderList, nei)
								end
							end
						end
					end
					createPCIR(point)
					for repCore in keys(clusters)
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
	diff = setdiff(1:num_records, rangeQueryPerformed)
	println("Unprocessed Points : ", length(diff), "\n")
	println("No. of Clusters : ", length(clusters), "\n")
end

function main(argv)
	println("\n-----AnyDBC Start------\n")
	s_time = time()
	file_name = argv
	global dataSet = readData(file_name)
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
