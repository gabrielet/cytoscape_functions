library("igraph")

#source("igraphManuale.R")

##################################################################
# using IGraph to mimic Centiscape and Pesca, two Cytoscape apps #
##################################################################

# this function allows to find nodes that belong
# to those shortest paths connecting two components
# input is a reference network and a list of edges
findNodesInPaths <- function(refNet, edgelist) {

	# declare empty vector
	nodes <- vector()
	# how many paths are in the edgelist?
	nPaths <- length(edgelist)
	# for each path
	for(i in 1:nPaths){
		# find nodes names
		nodes <- append(nodes, V(refNet)$name[unlist(edgelist[[i]])])
	}
	# return unique nodes in the shortest paths
	return(unique(nodes[-1]))
}

# find first neighbours for a set of nodes, i.e. the probe
# in the specified network
firstNeighbours <- function(network, probe) {
	# find first neighbours
	first <- induced.subgraph(network, unlist(neighborhood(network, order=1, probe)))
	return(first)
}

# this function is used to find all nodes in a network
# that are disconnected, i.e. do not have neighbours
# it is used after running `decompose` on a network
findNotConnectedNodes <- function(decomposednet) {
	# how many disconnected components are there?
	cmps <- length(decomposednet)
	# declare empty vector
	notConnected <- vector()

	# TODO FIND BIGGEST COMPONENT AND RETURN THE OTHERS!!!!!!!!!!!!!!!!

	# THIS IS NOT USING THE FIRST COMPONENT ASSUMING THAT'S THE BIGGEST ONE.
	# BUT IS THIS ALWAYS TRUE??
	for (i in 2:cmps) {
		notConnected <- c(notConnected, V(decomposednet[[i]])$name)
	}
	#senza la posizione uno!
	return(unique(notConnected[-1]))
}

#leggo il file sif
loadSif <- function(filename){
	sifFile <- read.table(filename)

	#creo la rete usando le colonne uno e tre del file sif
	tmpNet <- graph.data.frame(sifFile[,c(1,3)],directed=FALSE)

	#aggiusto la componente connessa, rimuovendo doppioni e loops
	network <- simplify(tmpNet, remove.multiple=TRUE, remove.loops=TRUE)

	if (length(decompose(network))>1) {
		print("verify that this network does not contain any disconnected component. if it does try spCluster or connectIsolated")
	}
	
	return(network)
}

#########
# PESCA #
#########

#probe contiene la lista dei nodi che voglio investigare: sono connessi o no?

spCluster <- function(referenceNet, probe) {
	#creo la sottorete e ne calcolo le componenti connesse e non
	subnet <- induced.subgraph(referenceNet, probe)
	decomposednet <- decompose(subnet)
	if (length(decomposednet)>1) {
		giant <- V(decomposednet[[1]])$name
		notCnctd <- findNotConnectedNodes(decomposednet)
		#sp-cluster
		sp <- get.all.shortest.paths(referenceNet, giant, notCnctd)
	
		#crea la lista di nodi contenuti negli short paths estratti
		nodeslist <- findNodesInPaths(referenceNet, sp$vpath)

		#aggiungo le due componenti, gigante e disconnessa
		nodeslist <- unique(c(nodeslist, giant, notCnctd))

		#e poi costruisco la rete, se uso la funzione corrente
		connectedNet <- induced.subgraph(referenceNet, nodeslist)
		return(connectedNet)
	}
	else{
		return("no disconnected component")
	}
}

connectIsolated <- function(referenceNet, probe) {
	#creo la sottorete e ne calcolo le componenti connesse e non
	subnet <- induced.subgraph(referenceNet, probe)
	decomposednet <- decompose(subnet)
	if (length(decomposednet)>1) {
		giant <- V(decomposednet[[1]])$name
		notCnctd <- findNotConnectedNodes(decomposednet)
	
		#connect disconnected component
		connect <- get.shortest.paths(referenceNet, giant, notCnctd)

		#crea la lista di nodi contenuti negli short paths estratti
		nodeslist <- findNodesInPaths(referenceNet, connect$vpath)

		#aggiungo le due componenti, gigante e disconnessa
		nodeslist <- unique(c(nodeslist, giant, notCnctd))

		#e poi costruisco la rete, se uso la funzione corrente
		connectedNet <- induced.subgraph(referenceNet, nodeslist)
		simplyNet <- simplify(connectedNet, remove.multiple=T, remove.loops=T)
		return(simplyNet)
	}
	else{
		return("no disconnected component")
	}
}

##############
# CENTISCAPE #
##############

#deg <- degree(connectedNet)
#ecc <- eccentricity(connectedNet)
#bet <- betweenness(connectedNet)
#clo <- closeness(connectedNet)

#alpha.centrality	Find Bonacich alpha centrality scores of network positions
#authority.score	Kleinberg's authority centrality scores.
#betweenness.estimate	Vertex and edge betweenness centrality
#centralization.evcent	Centralize a graph according to the eigenvector centrality of vertices
#centrEigen	Centralize a graph according to the eigenvector centrality of vertices
#closeness.estimate	Closeness centrality of vertices
#edge.betweenness	Vertex and edge betweenness centrality
#edge.betweenness.estimate	Vertex and edge betweenness centrality
#edgeBetweenness	Vertex and edge betweenness centrality
#eigenCentrality	Eigenvector Centrality Scores of Network Positions
#hub.score	Kleinberg's hub centrality scores.
#powerCentrality	Find Bonacich Power Centrality Scores of Network Positions
#subgraph.centrality	Find subgraph centrality scores of network positions

#network is the graph we want to weigh in order to perform the simulation
#simulations is the number of networks the function will generate
#[valMin-valMax] is the range in between the weight vectors will be generated
#centrality could be a single centrality or a vector of centralities (i.e. alphaCentrality,authorityScore,betweenness,closeness,edge.betweenness,eigenCentrality)

#the idea behind the algorithm: simulate a network with different topologies that depends from a weights vector that gives a multiplication factor to each node

#how it works:
#rgSim: for each simulation a weight vector is randomly created in between a range defined by the user. the number of simulation is also defined by the user. the input network is also defined by the user. it can be both directed and undirected
#addEdges: starting from the network and a random vector, the edgelist is visited: for each node the algorithm performs a loop for each copy the node has: the node 1 has 0 copies: 0 loops, the node 4 has 3 copies: 3 loops. For each loop a new node is added and a new set of edges is created. First the algorithm finds the position of the original node in the edge list. This step permits to get the neighbors of the original node that will be neighbors of each copy (rows: 33-42). Then an edge is added between the node and its copy (row: 43). Another edge is added if the network is considered as directed: the edges between the original node and its copies are not directed so a second edge A,B and B,A is added (rows: 44-46). Finally the new edges are added to the edgeList and the loops begin again with a new copy, if it exists.

addEdges <- function(network, weights, direction=FALSE){ #input: network, the weights and if the network is directed. default: undirected network

	allTheEdges <- get.edgelist(network) #the edgelist of the network
	newEdges <- c()
	len <- nrow(allTheEdges)
	newNode <- length(V(network)) #the new node i should add is nodes + 1
	copies <- weights - 1 #total number of copies each node will add to the network
	
	for(n in V(network)){ #for the node n
		if(copies[n]>0){ #if there is at least a copy
			for(cp in 1:copies[n]){ #for each copy				
				newNode <- newNode + 1 #a new node is going to be added
				pos <- which(allTheEdges==n) #find the node in the edgeList
				newEdges <- c()
				for(p in pos){	#for each position create a new couple of values: (n,oldNeighborOfNFather)
					if(p<=len){
						newEdges <- rbind(newEdges,c(newNode,allTheEdges[p+len]))
					}
					else{
						newEdges <- rbind(newEdges,c(allTheEdges[p-len],newNode))
					}
				}
				newEdges <- rbind(newEdges,c(n,newNode)) #add an edge between the new node and its father
				if(direction==TRUE){
					newEdges <- rbind(newEdges,c(newNode,n))
				}
				allTheEdges <- rbind(allTheEdges,newEdges)
				len <- nrow(allTheEdges) #how many edges
			}
		}
	}
	#COMMENT THESE TWO LINES
	newNetwork <- graph(c(t(allTheEdges)),directed=direction) #the new network
	return(newNetwork)
	#IF THE USERS NEEDS AN EDGE LIST INSTEAD OF A NETWORK THEN UNCOMMENT THE FOLLOWING AND COMMENT THE TWO LINES ABOVE
	#return(allTheEdges)
}

############################################################
#use this method when the network nodes has specific names!#
############################################################

addEdgesWithNames <- function(network, weights, direction=FALSE){ #input: network, the weights and if the network is directed. default: undirected network

	allTheEdges <- get.edgelist(network) #the edgelist of the network
	newEdges <- c()
	len <- nrow(allTheEdges)
	nnod <- length(V(network))
	newNode <- length(V(network)) #the new node i should add is nodes + 1
	copies <- weights - 1 #total number of copies each node will add to the network
	
	for(n in 1:nnod){ #for each original node
		if(copies[n]>0){ #if there is at least a copy
			for(cp in 1:copies[n]){ #for each copy				
				newNode <- paste(V(network)$name[n],cp,sep="_") #a new node is going to be added
				pos <- which(allTheEdges==V(network)$name[n]) #find the node in the edgeList
				newEdges <- c()
				for(p in pos){	#for each position create a new couple of values: (n,oldNeighborOfNFather)
					if(p<=len){
						newEdges <- rbind(newEdges,c(newNode,allTheEdges[p+len]))
					}
					else{
						newEdges <- rbind(newEdges,c(allTheEdges[p-len],newNode))
					}
				}
				newEdges <- rbind(newEdges,c(V(network)$name[n],newNode)) #add an edge between the new node and its father
				
				if(direction==TRUE){
					newEdges <- rbind(newEdges,c(newNode,n))
				}
				allTheEdges <- rbind(allTheEdges,newEdges)
				len <- nrow(allTheEdges) #how many edges
			}
		}
	}
	#COMMENT THESE TWO LINES
	newNetwork <- graph(c(t(allTheEdges)),directed=direction) #the new network
	return(newNetwork)
	#IF THE USERS NEEDS AN EDGE LIST INSTEAD OF A NETWORK THEN UNCOMMENT THE FOLLOWING AND COMMENT THE TWO LINES ABOVE
	#return(allTheEdges)
}


centralityPos <- function(nNodes,weight){ #finding the position that correspond to a node and its copies in the centrality vector
		
	copies <- weight-1
	maxval <- max(weight) #find the max multiplication factor (the maximal number of copies)
	whatSum <- rep(0,maxval)
	totalPos <- data.frame()
	pos <- 1
	for(i in copies){ #how many copies a node has?
		whatSum <- rep(0,maxval)
		if(i == 0){ #if the node has 0 copies the centrality is just itself
			whatSum[1] <- pos #in pos[1] will put the position of the original node in the centrality vector
			pos <- pos + 1
		}
		else{ #otherwise it should sum the centralities of all the copies and the original
			whatSum[1] <- pos #the pos[1] is always for the real position of the original node
			for(j in 1:copies[pos]){ #then for the other copies
				nNodes <- nNodes + 1 #copy j
				whatSum[j+1] <- nNodes #add the position of the copy j to the vector of this original node
			}
			pos <- pos + 1 #go for the next node
		}
		totalPos <- rbind(totalPos,whatSum) #add a row, one for each node
	}
	return(totalPos) #the final matrix has nrow as the node in the original network and ncol has the maximum number of copies a node has in the network
}	

centralitySum <- function(centralityTmp,weight){ #sum up the centralities (compute the centralities of the multiplied network in "centralityTmp") of each original node and its copies
	
	centrality <- rep(0,length(weight)) #the centrality vector has the same lenght of the number of orignal nodes
	tmp <- centralityPos(length(weight),weight) #computing the position of the centralities to sum for each original node
	
	for(i in 1:nrow(tmp)){ #for each row (a node)
		for(j in 1:ncol(tmp[1,])){ #for each column (the original node and its copies)
			pos <- unlist(tmp[i,][j]) #get the position of the centrality to add to total sum
			if(length(centralityTmp[pos])!=0){ #if there actually is a centrality to add (there is at least a copy)
				centrality[i] <- centrality[i] + centralityTmp[pos] #sum the correspondend value
			}
		}
	}
	return(centrality)
}

rgSim <- function(network, simulations, valMin, valMax, centrality="closeness", direction=FALSE){

	weight <- vector()
	centralities <- vector()
	clo <- vector()
	betw <- vector()
	alpha <- vector()
	vertexes <- vcount(network)
	meansClo <- vector()
	meansAlpha <- vector()
	meansBetw <- vector()	
	label <- c("closeness","alpha.centrality","betweenness")
	results <- data.frame()
		
	if(all(centrality %in% label, na.rm=FALSE)){#if all the centralities are correctly typed then go on
		for(i in 1:simulations){#for each simulation an vector of random weights is generated
			weight <- round(runif(vertexes,valMin,valMax)) #generates #vertexes weights within valMax and valMin
			#add the new edges to the network and then compute the centralities
			
			#IF THE ADDEDGES FUNCTION RETURNS AN EDGE LIST THEN UNCOMMENT THESE TWO LINES ...
			#newEdgesList <- addEdges(network,weight)			
			#improvedNetwork <- graph(c(t(newEdgesList)))
			#AND COMMENT THE LINE BELOW			
			improvedNetwork <- addEdges(network,weight,direction)
			#plot(improvedNetwork)dev.new() #to print the multiplied networks, uncomment the line.
			
			if("closeness" %in% centrality){
				centrTmp <- closeness(improvedNetwork)
				centralitiesClo <- centralitySum(centrTmp,weight)
				clo <- c(clo,centralitiesClo)
				#meansClo[i] <- mean(centralitiesClo) #compute the mean closeness for the simulation #i
			}
			if("alpha.centrality" %in% centrality){
				centrTmp <- alpha.centrality(improvedNetwork)
				centralitiesAlpha <- centralitySum(centrTmp,weight)
				alpha <- c(alpha,centralitiesAlpha)
				#meansAlpha[i] <- mean(centralitiesAlpha)
			}
			if("betweenness" %in% centrality){
				centrTmp <- betweenness(improvedNetwork)
				centralitiesBetw <- centralitySum(centrTmp,weight)
				betw <- c(betw,centralitiesBetw)
				#meansBetw[i] <- mean(centralitiesBetw)
			}
		}
	}
	else{
		stop("One of the centralities you selected does not exist")
	}
	#results <- c(mean(meansClo),mean(meansBetw),mean(meansAlpha))
	row <- (length(V(network))*simulations)
	results <- matrix(nrow=row,ncol=3)
	if("closeness" %in% centrality){results[,1] <- t(clo[-1])}else{results[,1] <- t(vector(0,row))}
	if("alpha.centrality" %in% centrality){results[,2] <- t(alpha[-1])}else{results[,2] <- t(vector(0,row))}
	if("betweenness" %in% centrality){results[,3] <- t(betw[-1])}else{results[,3] <- t(vector(0,row))}
	colnames(results) <- label
	return(results)
}

createWeights <- function(howmany,min,max){
	
	valMin <- round(min)
	valMax <- round(max)
	weight <- round(runif(howmany,valMin,valMax)) #generates howmany weights within valMax and valMin
	return(weight)
}

#some utilities for network analysis

#####################################################
#aCsv: must be a csv table which contains information about the nodes in the graph and their eccentricity
#anotherCsv: must be a csv table which contains information about the nodes in the graph and their eccentricity
#centrality: is a string which represents the name of the attribute relative to eccentricity
#flag: does the csv have an header?

#return a list of nodes, as defined by the IGraph library, which shows a difference in the Eccentricity centrality between the two graphs
#####################################################

#this function allows the user to use as input two loaded csv files
eccentricityMiningFromCsv <- function(aGraph, anotherGraph, centrality, flag=TRUE){
	if(centrality %in% colnames(aGraph) && centrality %in% colnames(anotherGraph)){
		diffNodes <- vector()
		for(i in aGraph$name){
			aPos <- which(aGraph$name==i)#get the position of the node i in aGraph
			anotherPos <- which(anotherGraph$name==i)#the same for anotherGraph
			aTmp <- aGraph[,centrality][aPos]#get its eccentricity
			anotherTmp <- anotherGraph[,centrality][anotherPos]#again, for anotherGraph
			if(length(aTmp)){#if aTmp contains something
				if(length(anotherTmp)){#and anotherTmp too
					if(round(aTmp,4)!=round(anotherTmp,4)){#then store the name of the node
						diffNodes <- append(i,diffNodes)
					}
				}
			}
		}
		return(diffNodes)
	}
	else{
		print("chosen centrality is not an existing column. Check these lists:")
		print(colnames(aGraph))
		print(colnames(anotherGraph))
	}
}

#this function allows the user to use as input two strings which represents the filenames
eccentricityMiningFromFilename <- function(aCsv, anotherCsv, centrality, flag=TRUE){

	aGraph <- read.csv(aCsv,header=flag)
	anotherGraph <- read.csv(anotherCsv,header=flag)
	if(centrality %in% colnames(aGraph) && centrality %in% colnames(anotherGraph)){
		diffNodes <- vector()
		if(file.exists(aCsv) && file.exists(anotherCsv)){
			for(i in aGraph$name){
				aPos <- which(aGraph$name==i)#get the position of the node i in aGraph
				anotherPos <- which(anotherGraph$name==i)#the same for anotherGraph
				aTmp <- aGraph[,centrality][aPos]#get its eccentricity
				anotherTmp <- anotherGraph[,centrality][anotherPos]#again, for anotherGraph
				if(length(aTmp)){#if aTmp contains something
					if(length(anotherTmp)){#and anotherTmp too
						if(round(aTmp,4)!=round(anotherTmp,4)){#then store the name of the node
							diffNodes <- append(i,diffNodes)
						}
					}
				}
			}
			return(diffNodes)
		}
		else print("one of the files does not exist")
	}
	else{
		print("chosen centrality is not an existing column. Check these lists:")
		print(colnames(aGraph))
		print(colnames(anotherGraph))
	}
}


#####################################################
#aGraph is a graph as defined by the IGraph library
#aNode is a node as defined by the IGraph library

#return a list of paths
#####################################################
findLongestPath <- function(aGraph,aNode){
	
	paths <- get.all.shortest.paths(aGraph,from=aNode,to=V(aGraph),mode="all")
	tmp <- 1
	intPaths <- list()
	#find the length of the longest shortest path generated from aNode
	for(i in 1:length(paths$res)){
		l <- length(paths$res[[i]])
		if(l>tmp){
			tmp <- l
		}
	}
	#now check if there are more paths which have that length
	j <- 1
	for(i in 1:length(paths$res)){
		if(length(paths$res[[i]])==tmp){
			intPaths[[j]] <- paths$res[[i]]
			j <- j+1
		}
	}
	return(intPaths)
}


#####################################################
#targets is a file which contains the proteins to search
#paths is a list of paths

#return a list of paths associated to the targets
#####################################################
findTargets <- function(targets, paths){
	
	l <- 1
	targetedPaths <- list()
	if(file.exists(targets)){
		file <- read.table(targets)
		for(i in 1:length(aPaths)){#how many lists are in paths?
			for(j in 1:length(aPaths[[i]])){#for each list, how many objects are stored?
				for(k in 1:dim(file)[1]){#search a protein in all the objects
					if(file[k,] %in% aPaths[[i]][[j]]$name){
						targetedPaths[[l]] <- append(targetedPaths,toString(aPaths[[i]][[j]]$name))
						l <- l+1
					}
				}
			}
		}
	}
	else{
		print("the chosen file does not exist")
	}
	return(targetedPaths)
}

#the same but with the specification doesNotContain a set of proteins
findTargetsWithoutSet <- function(contains, notContains, paths){
	
	l <- 1
	targetedPaths <- list()
	if(file.exists(targets)){
		contains <- read.table(targets)
		
		for(i in 1:length(aPaths)){#how many lists are in paths?
			for(j in 1:length(aPaths[[i]])){#for each list, how many objects are stored?
				for(k in 1:dim(contains)[1]){#search a protein in all the objects
					if(file[k,] %in% aPaths[[i]][[j]]$name){
						targetedPaths[[l]] <- append(targetedPaths,toString(aPaths[[i]][[j]]$name))
						l <- l+1
					}
				}
			}
		}
	}
	else{
		print("the chosen file does not exist")
	}
	return(targetedPaths)
}
