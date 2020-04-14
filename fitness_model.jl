# Import required packages
import Pkg
import Conda
conda.add("matplotlib"; channel="conda-forge")

Pkg.add("BioSequences")
Pkg.add("BlockArrays")
Pkg.add("Colors")
Pkg.add("GraphPlot")
Pkg.add("LightGraphs")
Pkg.add("LinearAlgebra")
Pkg.add("Plots")
Pkg.add("PyCall")
Pkg.add("Statistics")
Pkg.add("StatsBase")
Pkg.add("StatsPlots")

using BioSequences
using BlockArrays
using Colors
using GraphPlot
using LightGraphs
using LinearAlgebra
using Pkg
using Plots
pyplot()
using Statistics
using StatsBase
using StatsPlots




#################################################
# Make Dictionary or Tuples list of Nodes
#################################################

"""
makeNodeTuples(x)

Creates a list of tuples made up of all unique combinations of distinct nucleotides
in a 4ᴸ RNA chain.  

"""
function makeNodeTuples(L)
    nucleotides = ['A','U','C','G'] # list of possible RNA nucleotides
    nodes = reverse.(Iterators.product(fill(nucleotides,L)...))[:]; #Create tuple of all possible combinations
    return nodes
end

"""
makeNode(k,x)

Creates a string by taking x tuple's initial nucleotide charachter and concatenates 
all the subsequent characters in the tuple - this function is only used by the function
makeNodeList.  

"""
function makeNode(nodes, L,x)
    node = nodes[x][1]
    for i = 2:L
        node = node*nodes[x][i]
    end
    return node
end

"""
makeNodeList(nodes)

Creates an array of all the possible nucleotide combinations.  

"""
function makeNodeList(nodes)
    L = convert(Int,log(4,length(nodes)))
    nodeList=String[]
    for i in eachindex(nodes)
        x = makeNode(nodes,L,i)
           push!(nodeList,x)
    end
    return nodeList
end

"""
addFitnessValues(fileName)

import file of fitness values return a dictionary of nucleotide and fitness value.  

"""
function addFitnessValues(fileName)
    io = read(open(fileName),String);
    io = replace(io[4:end], "\r\n"=>";");
    io = replace(io[1:end], "T"=>"U");
    
    fvDict = Dict()
    fvArray = split(io,";")
    
    for i in eachindex(fvArray)
        x = split(fvArray[i],",")
        tempDict=Dict(x[1] => parse(Float64, x[2]))
        fvDict = merge(fvDict, tempDict)
    end
    fvDict = sort(collect(fvDict),by=x->x[1])
    return fvDict
end

function graphExampleCircle(L)
    nodes₁ = makeNodeTuples(L);
    adjMatrix₁,null = makeAdjMatrix(nodes₁)
    graph₁=Graph(adjMatrix₁)
    gplot(graph₁,nodefillc="orange", layout=circular_layout, nodelabel=makeNodeList(nodes₁))
end

function graphExample(L)
    nodes₁ = makeNodeTuples(L);
    adjMatrix₁,null = makeAdjMatrix(nodes₁)
    graph₁=Graph(adjMatrix₁)
    gplot(graph₁,nodefillc="orange", nodelabel=makeNodeList(nodes₁))
end

#################################################
# Make Adjacency Matrices
#################################################

"""
makeAdjMatrix(nodes)

Creates an adjacency matrix of all the possible nucleotide combinations that are one 
mutation away from each other.  

"""
function makeAdjMatrix(nodes)
    
    L = convert(Int,log(4,length(nodes)))
    
    adjMatrix = zeros(Int64,4^L,4^L)
    rndMatrix = zeros(Int64,4^L,4^L)
    n = size(adjMatrix,1)

    for i in 1:n
        for j in 1:n
            if count(nodes[i].!= nodes[j]) == 1
                adjMatrix[i,j] = 1 
            end
            if count(nodes[i].!= nodes[j]) == 2
                rndMatrix[i,j] = 1 
            end   
        end
    end
    return adjMatrix, rndMatrix
end

"""
makeFitAdjMatrix(nodes)

Creates an adjacency matrix of all the possible nucleotide combinations that are one 
mutation away from each other.  

"""
function makeFitAdjMatrix(nodes)
    
    L = convert(Int,log(4,length(nodes)))
    fitMatrix = zeros(Float64,4^L,4^L)
    adjMatrix = zeros(Float64,4^L,4^L)
    n = size(adjMatrix,1)
    
    for i in 1:n
        for j in 1:n
            
            counter = 0
            
            for k in 1:L
                if  nodes[i][1][k] != nodes[j][1][k]
                    counter = counter + 1
                end
            end
            
            if counter >= 2 || counter == 0
                adjMatrix[i,j] =  0
                fitMatrix[i,j] =  0
            else
                adjMatrix[i,j] =  1
                fitMatrix[i,j] =  1 + nodes[i][2] - nodes[j][2]
            end
            
        end
    end
    combinedMatrix = adjMatrix.*fitMatrix
    return combinedMatrix
end

"""
makeRandomMutation(nodes)

Add a random mutation that are two mutations away from each other.  

"""
function makeRandomMutation(rndMatrix)
    
    L = convert(Int,log(4,size(rndMatrix,1)))
    
    randRange = findall(!isequal(0),rndMatrix[:,:])
    jump = rand(randRange)
    
    addRndMatrix = zeros(Int64,4^L,4^L)
    addRndMatrix[jump[1],jump[2]] = 1
    addRndMatrix[jump[2],jump[1]] = 1
    
    return addRndMatrix
end

#################################################
# Make Paths
#################################################

"""
makePath(startNode, endNode)

Creates an array of a random path between the startNode and endNode

"""
function makePath(startNode, endNode, nodeList, adjMatrix)
    pathList = String[];
    jump = startNode;
    push!(pathList,nodeList[jump])

   while jump != endNode
        randRange= findall(!isequal(0),adjMatrix[jump,:]) # find the locations of all the possible paths from the current node
        jump = rand(randRange) # randomly select from all the possible paths to the next node
        push!(pathList,nodeList[jump]) # push the current node we're evaluating from on to the path list.
    end
    return pathList
end

"""
identifyPaths(iter, startNode, endNode, nodeList, rndMatrix,adjMatrix, rndEdge)

Creates arrays of counts and paths for a 'iter' number of random path between 
the startNode and endNode

"""
function identifyPaths(iter, startNode, endNode, nodeList, rndMatrix, adjMatrix, randEdge = 0)

    pathCount =  Array{Int}(undef,iter)
    pathList =  Array{String}(undef,iter)

    if randEdge == 1
        for i in 1:iter
            addRndMatrix = makeRandomMutation(rndMatrix)
            path = makePath(startNode,endNode,nodeList,adjMatrix+addRndMatrix)
            pathCount[i] = length(path)
            pathList[i] = join(path) 
        end    
    elseif randEdge == 0
        for i in 1:iter
            path = makePath(startNode,endNode,nodeList, adjMatrix)
            pathCount[i] = length(path)
            pathList[i] = join(path)
        end
    end

    pathMedian = median(pathCount)
    pathMax = maximum(pathCount)
    pathMin = minimum(pathCount)
    println("Median Path Length: ",pathMedian)
    println("Maximum Path Length: ",pathMax)
    println("Minimum Path Length: ",pathMin)

    return pathCount, pathList
end

"""
makeFitPath(startNode, endNode, nodeList, adjMatrix)

Creates an array of a random path between the startNode and endNode

"""
function makeFitPath(startNode, endNode, nodeList, adjMatrix,level)
    pathList = String[];
    jump = startNode;
    push!(pathList,String(nodeList[jump][1]))

   while jump != endNode 
        randRange= findall(x -> x >= level,adjMatrix[jump,:]) # find the locations of all the possible paths from the current node
        
        if isempty(randRange)
            push!(pathList,String("XXXXXXX"))
            break
        end
        
        jump = rand(randRange) # randomly select from all the possible paths to the next node
        push!(pathList,String(nodeList[jump][1])) # push the current node we're evaluating from on to the path list.
    end
    return pathList
end

"""
identifyFitnessPaths(iter, startNode, endNode, nodeList, fitMatrix)

Creates arrays of counts and paths for a 'iter' number of random path between 
the startNode and endNode

"""
function identifyFitnessPaths(iter, startNode, endNode, nodeList, adjMatrix, level)

    pathCount =  Array{Int}(undef,iter)
    pathList =  Array{String}(undef,iter)

    for i in 1:iter
        path = makeFitPath(startNode,endNode,nodeList,adjMatrix,level)
        pathCount[i] = length(path)
        pathList[i] = join(path)
    end

    println("Median Path Length: ", median(pathCount))
    println("Maximum Path Length: ",maximum(pathCount))
    println("Minimum Path Length: ",minimum(pathCount))

    return pathCount, pathList
end

#################################################
# Identify Successful Nodes
#################################################
"""
findShortPathNodes(pathList, pathCount)

Create a dictionary of all the nodes on all successful paths

"""
function findShortPathNodes(pathList,pathCount,L)
    # Find top 10% shortest paths
    orderedPathList = sort(pathList,by=length);
    topOrderedPathList = orderedPathList[1:Int(length(pathCount)/10)];

    # Find top 10% most frequent path size
    freqPathCount = countmap(pathCount);
    sortFreqPathCount = sort(collect(freqPathCount), by=x->x[2], rev = true);
    topOrderedPathCount = sortFreqPathCount[1:round(Int,size(sortFreqPathCount)[1]/10)];
    
    # Find the size of the Array to build
    counter1=0;
    for i in eachindex(topOrderedPathCount)
        counter1 = counter1 +topOrderedPathCount[i][2]
    end
    
    # Find the paths of the 10% most frequent size
    topPathSizeSeq = Array{String}(undef,counter1)
    counter2 = 1
    for i in eachindex(topOrderedPathCount)
        for j in eachindex(pathList) 
            if length(pathList[j])== (L*topOrderedPathCount[i][1])
               topPathSizeSeq[counter2] = pathList[j]
               counter2 = counter2+1 
            end
        end
    end
    
    sumTotal = 0;
    for i in eachindex(topPathSizeSeq)
            sumTotal = sumTotal + length(topPathSizeSeq[i])
    end
    
    list = Array{String}(undef,Int(sumTotal/L))
    counter3 = 1
    for i in eachindex(topPathSizeSeq)
        for j = 1:L:length(topPathSizeSeq[i])-L+1
            list[counter3] = join(topPathSizeSeq[i][j:j+(L-1)])
            counter3 = counter3 + 1
        end
    end
    sortedTopSeqeunces = sort(collect(countmap(list)), by=x->x[2], rev=true)

    nodeCountlist = Array{Int}(undef,length(sortedTopSeqeunces))
    for i in eachindex(sortedTopSeqeunces)
        value = sortedTopSeqeunces[i][2]
        nodeCountlist[i] = value
    end

    return sortedTopSeqeunces, nodeCountlist
end

"""
shortPath(pathList, pathCount,L)

Create an array of sequences on the shortest path
"""
function shortPath(pathList, pathCount, L)
    pathString = pathList[findall(x -> x == minimum(pathCount),pathCount)]
    shortPath = Array{String}(undef, Int(length(pathString[1])/L))
    jj=0
    for i in eachindex(shortPath)
        shortPath[i] = pathString[1][1+(L*jj):7+(L*jj)]
        jj +=1
    end
    return shortPath
end