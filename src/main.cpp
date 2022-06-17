#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <time.h>
#include "../include/match.hpp"
#include "../include/argloader.hpp"
#include "../include/argraph.hpp"
#include "../include/argedit.hpp"
#include "../include/nodesorter.hpp"
#include "../include/probability_strategy.hpp"
#include "../include/nodesorter.hpp"
#include "../include/nodeclassifier.hpp"

#include "../include/vf3_sub_state.hpp"

#include <vector>
using namespace std;


typedef string data_t;                                                  // node attributes are strings
#define PARAMETER_K_IF_NEEDED
template<> long long VF3SubState<data_t,data_t,Empty,Empty>::instance_count = 0;

typedef vector<pair<data_t, data_t>> MU;
typedef vector<MU> SOL;
typedef vector<pair<unsigned int, SOL>> SOLUTION;

typedef struct visitor_data_s
{                                   
  MU mu;
  SOL sol;
} visitor_data_t;


struct TreeNode
{
	TreeNode *parent;
	vector<TreeNode *> children;
	MU value;
	unsigned int pg_index;
	unsigned int dg_index;
};

vector<ARGraph<data_t, Empty>> G, H;                                  // vectors containing snapshots for G and H


bool visitor(int n, node_id nH[], node_id nG[], void *state, void *usr_data)
{
	visitor_data_t* data = (visitor_data_t*)usr_data;

	VF3SubState<data_t, data_t, Empty, Empty>* s = static_cast<VF3SubState<data_t, data_t, Empty, Empty>*>(state);

	for(int k = 0; k < n; k++)                     // n: no of nodes in pattern graph
	{                                                     
		if(nH[k] != NULL_NODE) 
		{
			data->mu.push_back(make_pair(s->GetGraph1()->GetNodeAttr(nH[k]), s->GetGraph2()->GetNodeAttr(nG[nH[k]])));                         // complete solution
		}
	}

	data->sol.push_back(data->mu);
	data->mu.clear();

	return false;
}



// original vf3 for static graphs
bool VF3(ARGraph<data_t, Empty> H, ARGraph<data_t, Empty> G, SOL &soln)
{
	int n=0;
	visitor_data_t vis_data;

	int nodesH, nodesG;
	nodesH = H.NodeCount();
	nodesG = G.NodeCount();
	node_id *nH = new node_id[nodesH];
	node_id *nG = new node_id[nodesG];

	NodeClassifier<data_t, Empty> classifier(&G);
	NodeClassifier<data_t, Empty> classifier2(&H, classifier);

	vector<int> class_H = classifier2.GetClasses();
	vector<int> class_G = classifier.GetClasses();

	VF3NodeSorter<data_t, Empty, SubIsoNodeProbability<data_t, Empty> > sorter(&G);

	vector<node_id> sorted = sorter.SortNodes(&H);

	VF3SubState<data_t, data_t, Empty, Empty> s0(&H, &G, class_H.data(), class_G.data(), classifier.CountClasses(), sorted.data() PARAMETER_K_IF_NEEDED);
	match<VF3SubState<data_t, data_t, Empty, Empty> >(s0, &n, nH, nG, visitor, &vis_data);

	if (!vis_data.sol.empty())
	{
    		soln = vis_data.sol;
    		return true;
    	}
    	else
    		return false;
}





void MapFirstPatternGraph(unsigned int i, SOLUTION &solutioni, TreeNode *root, int window_size_increment) 
{
	double timeAllMFPG = 0;
	unsigned long ticksMFPG = 0;
	unsigned long ticksFinalMFPG = 0;

	ticksMFPG = clock();

	SOL solutionij;

	for (unsigned int j = 0; j <= G.size() - H.size() + window_size_increment && j < G.size(); j++)         // dont apply first VF3 for every target graph snapshot
	{                  
		if (VF3(H[i], G[j], solutionij) == true) 
		{
			solutioni.push_back(make_pair(j,solutionij));

			for(SOL::iterator itsol = solutionij.begin(); itsol != solutionij.end(); ++itsol) 
			{

				// make new tree node for every mu
				TreeNode *node = new TreeNode;

				//give value to tree node
				node->parent = root;
				node->value = *itsol;                  // node contains complete mu
				node->pg_index = i;
				node->dg_index = j;
				root->children.push_back(node);
				
			}
		}
	}

	ticksFinalMFPG = clock() - ticksMFPG;
	timeAllMFPG = (double) ticksFinalMFPG/CLOCKS_PER_SEC;
	//cout << "time for mapping all graph pairs in MapFirstPatternGraph:  " << timeAllMFPG << endl;

}


// new vf3 for dynamic graphs
bool VF3(ARGraph<data_t, Empty> Hi_1, ARGraph<data_t, Empty> Gj, ARGraph<data_t, Empty> Hi, ARGraph<data_t, Empty> Gf, MU mu, vector<node_id> F_Hi, vector<node_id> F_Gf, SOL &soln)
{

	// mu corresponds to graphs Hi_1, Gj

	int nodesHi_1, nodesGj, nodesHi, nodesGf;
	nodesHi_1 = Hi_1.NodeCount();
	nodesGj = Gj.NodeCount();
	nodesHi = Hi.NodeCount();
	nodesGf = Gf.NodeCount();
	node_id *nH = new node_id[nodesHi];
	node_id *nG = new node_id[nodesGf];

	vector<node_id> nodesinmuHi, nodesinmuGf;

	int n=0;
	visitor_data_t vis_data;

	for(int p = 0; p < mu.size(); p++)
	{
		if(p < nodesHi && p < nodesHi_1)
		{
			nodesinmuHi.push_back(Hi.GetNodeId(mu[p].first));
		}
		if(p < nodesGj && p < nodesGf)
		{
			nodesinmuGf.push_back(Gf.GetNodeId(mu[p].second));
		}
	}

	NodeClassifier<data_t, Empty> classifier(&Gf);
	NodeClassifier<data_t, Empty> classifier2(&Hi, classifier);

	vector<int> class_Hi = classifier2.GetClasses();
	vector<int> class_Gf = classifier.GetClasses();

	VF3NodeSorter<data_t, Empty, SubIsoNodeProbability<data_t, Empty> > sorter(&Gf);

	vector<node_id> sorted = sorter.SortNodes(&Hi);

	vector<node_id> sorted1;
	for(auto it = sorted.begin(); it != sorted.end(); ++it)
	{
		auto iter = find(nodesinmuHi.begin(), nodesinmuHi.end(),*it);
		if(iter == nodesinmuHi.end())
		{
			sorted1.push_back(*it);
		}
	}

	// zero padding for sorted
	for(unsigned int s = 0; s < nodesinmuHi.size(); s++)
	{
		sorted1.insert(sorted1.begin(), 0);
	}


	// make initial state using nodesinmuHi and nodesinmuGf
	VF3SubState<data_t, data_t, Empty, Empty> s(&Hi, &Gf, class_Hi.data(), class_Gf.data(), classifier.CountClasses(), sorted1.data(), nodesinmuHi, nodesinmuGf, F_Hi, F_Gf);
	match<VF3SubState<data_t, data_t, Empty, Empty> >(s, &n, nH, nG, visitor, &vis_data);

    	if (!vis_data.sol.empty())
    	{
    		soln = vis_data.sol;
    		return true;
    	}

    	else
    		return false;
}





void FindNewMapping(unsigned int i, vector<node_id> NAM, vector<node_id> NNM, SOLUTION &solutioni, SOLUTION &solutioni_1, vector<TreeNode *> &Li, vector<TreeNode *> &Li_1, int window_size_increment)
{
	
	double timeAllFNM = 0;
	unsigned long ticksFNM = 0;
	unsigned long ticksFinalFNM = 0;

	ticksFNM = clock();


	// feasibility set of Hi
	vector<node_id> F_Hi, F_Gf;
	SOL solutionif;


	/*
	 * Calculate F_Hi:
	 * both OutEdge and InEdge are used to calculate F_Hi
	 * for both DIRECTED and UNDIRECTED graphs
	 *
	 */
	for(auto itNNM = NNM.begin(); itNNM != NNM.end(); itNNM++) 
	{
		for(int k = 0; k < H[i].NodeCount(); k++) 
		{
			if(k < H[i].OutEdgeCount(*itNNM)) 
			{
				node_id n_out = H[i].GetOutEdge(*itNNM, k);
				if((find(NAM.begin(), NAM.end(), n_out) != NAM.end()) && (find(F_Hi.begin(), F_Hi.end(), *itNNM) == F_Hi.end()))
					F_Hi.push_back(*itNNM);
			}

			if(k < H[i].InEdgeCount(*itNNM))
			{
				node_id n_in = H[i].GetInEdge(*itNNM, k);
				if((find(NAM.begin(), NAM.end(), n_in) != NAM.end()) && (find(F_Hi.begin(), F_Hi.end(), *itNNM) == F_Hi.end()))
					F_Hi.push_back(*itNNM);
			}
		}
	}

	vector<data_t> nodesHi;
	for(node_id nHi = 0; nHi < H[i].NodeCount(); nHi++)
	{
		nodesHi.push_back(H[i].GetNodeAttr(nHi));
	}

	for(unsigned int j = i - 1; j <= G.size() - H.size() + i - 1 + window_size_increment && j < G.size(); j++)  // search for previously found solutions within the window corresponding to H_{i-1}
	{			
		for(unsigned int k = 0; k < solutioni_1.size(); k++) 
		{
			if(j == solutioni_1[k].first) 
			{
				
				SOL solu = solutioni_1[k].second;                      // contains all previous mappings     

				for(unsigned int a = 0; a < solu.size(); a++)
				{

					MU mue = solu[a];

					// take nodes present in H[i]; filter using node attributes
					MU mue_new;
					for(unsigned int it = 0; it < mue.size(); it++)
					{
						auto iter = find(nodesHi.begin(), nodesHi.end(), mue[it].first);
						if(iter != nodesHi.end())
						{
							mue_new.push_back(make_pair(*iter, mue[it].second));
						}
					}


					//cout << "mue_new:" << endl;           // mapping corresponding to G[j]; 
                                                                                // it's validity for G[f] is verified by finding the nodes
                                                                                // present in mue_new as well as all edges between them in G[f] 

					vector<data_t> mue_new_first, mue_new_second;
					for(unsigned int i = 0; i < mue_new.size(); i++)
					{	
						mue_new_first.push_back(mue_new[i].first);
						mue_new_second.push_back(mue_new[i].second);

					}	

					// check 
					// the 
					// future
					for(unsigned int f = j+1; f <= G.size() - H.size() + i + window_size_increment && f < G.size(); f++)   // sliding window to look into the future
					{					
						
						// nodes of graph G[f]
						vector<data_t> nodesGf;
						for(node_id nGf = 0; nGf < G[f].NodeCount(); nGf++)
						{
							nodesGf.push_back(G[f].GetNodeAttr(nGf));
						}


						// look for nodes present in mue_new in graph G[f]
						vector<pair<int, data_t>> nodesfound;
						vector<data_t> nodesfound_second;
						int index;
						for(unsigned int l = 0; l < mue_new.size(); l++) 
						{
							auto it = find(nodesGf.begin(), nodesGf.end(), mue_new[l].second);
							if(it != nodesGf.end()) {
								index = distance(nodesGf.begin(), it);
								nodesfound.push_back(make_pair(index, mue_new[l].second));
								nodesfound_second.push_back(mue_new[l].second);
							}
						}

						
						// look for edges between nodes present in mue_new in graph G[f]
						vector<vector<data_t> > neighsHi, neighsGf;
						vector<data_t> neighHi, neighGf;
						int mue_new_edgeCount_Hi = 0;
						int mue_new_edgeCount_Gf = 0;
						for(unsigned int k = 0; k < mue_new.size(); k++)
						{
							for(unsigned int l = 0; l < mue_new.size(); l++)
							{
								if(H[i].HasEdge(H[i].GetNodeId(mue_new[k].first), H[i].GetNodeId(mue_new[l].first)))
								{
									neighHi.push_back(mue_new[l].first);        // find the neighbour of the node of Hi in mue_new
								}

								if(G[f].HasEdge(G[f].GetNodeId(mue_new[k].second), G[f].GetNodeId(mue_new[l].second)))
								{
									neighGf.push_back(mue_new[l].second);       // find the neighbour of the node of Gf in mue_new
								}

							}

							mue_new_edgeCount_Hi += neighHi.size();    // number of edges between nodes present in mue_new in H[i] 
							mue_new_edgeCount_Gf += neighGf.size();    // number of edges between nodes present in mue_new in G[f] 

							if(!neighHi.empty())
								neighsHi.push_back(neighHi);
							if(!neighGf.empty())
								neighsGf.push_back(neighGf);

							neighHi.clear();
							neighGf.clear();

						}

						
						int index_neigh_Hi, index_neigh_Gf;
						vector<pair<data_t, data_t> > HiEdgesFound;

						for(unsigned int k = 0; k < neighsHi.size() && k < neighsGf.size(); k++)
						{
							for(unsigned int l = 0; l < neighsHi[k].size() && l < neighsGf[k].size(); l++)
							{								

								auto iter_Hi = find(mue_new_first.begin(), mue_new_first.end(), neighsHi[k][l]);
								if(iter_Hi != mue_new_first.end())
								{
									index_neigh_Hi = distance(mue_new_first.begin(), iter_Hi);	
								}

								auto iter_Gf = find(mue_new_second.begin(), mue_new_second.end(), neighsGf[k][l]);
								if(iter_Gf != mue_new_second.end())
								{
									index_neigh_Gf = distance(mue_new_second.begin(), iter_Gf);	
								}

								if(index_neigh_Hi == index_neigh_Gf)
								{
									HiEdgesFound.push_back(make_pair(mue_new_first[k], neighsHi[k][l]));	
								}

							}
						}



						/* 
						 * Calculate F_Gf:
						 * both OutEdge and InEdge are used to calculate F_Gf
						 * for both DIRECTED and UNDIRECTED graphs
						 *
						 */


						node_id n_ouut, n_iin, n1; 
						if(nodesfound.size() == mue_new.size() && nodesfound_second.size() == mue_new.size() && HiEdgesFound.size() == mue_new_edgeCount_Hi)   
						{   
							// all nodes in mue_new as well as all edges between those nodes are found in G[f]
						            
							for(int it = 0; it < nodesfound.size(); it++)
							{
								for(n1 = 0; n1 < G[f].NodeCount(); n1++)
								{
									if(n1 < G[f].OutEdgeCount(nodesfound[it].first))
									{
										n_ouut = G[f].GetOutEdge(nodesfound[it].first, n1);
										auto itnodesfound_second = find(nodesfound_second.begin(), nodesfound_second.end(), G[f].GetNodeAttr(n_ouut));
										auto itF_Gf = find(F_Gf.begin(), F_Gf.end(), n_ouut);
										if(itnodesfound_second == nodesfound_second.end() && itF_Gf == F_Gf.end())
											F_Gf.push_back(n_ouut);
									}

									if(n1 < G[f].InEdgeCount(nodesfound[it].first))
									{
										n_iin = G[f].GetInEdge(nodesfound[it].first, n1);
										auto itnodesfound_second = find(nodesfound_second.begin(), nodesfound_second.end(), G[f].GetNodeAttr(n_iin));
										auto itF_Gf = find(F_Gf.begin(), F_Gf.end(), n_iin);
										if(itnodesfound_second == nodesfound_second.end() && itF_Gf == F_Gf.end())
											F_Gf.push_back(n_iin);
									}

								}
							}
						}



						// Apply VF3 with non-empty initial state
						if(!F_Gf.empty()) 
						{
							if(VF3(H[i-1], G[j], H[i], G[f], mue_new, F_Hi, F_Gf, solutionif) == true)
							{
								solutioni.push_back(make_pair(f,solutionif)); 

								bool contains = false;
								for(SOL::iterator itsol=solutionif.begin(); itsol!=solutionif.end();++itsol)
								{
								
									TreeNode *node = new TreeNode;

									for(auto it = Li_1.begin(); it != Li_1.end(); it++)
									{
										if((*it)->value == mue && (*it)->dg_index == j)                     // to identify the correct parent node
										{                     
											node->parent = *it;                                         // find new mapping
											node->value = *itsol;
											node->pg_index = i;
											node->dg_index = f;

											for(auto iter = (*it)->children.begin(); iter != (*it)->children.end(); ++iter)
											{
												if((*iter)->value == *itsol && (*iter)->dg_index == f)
												{
													contains = true;
												}

											}
											if(!contains)
											{
												(*it)->children.push_back(node);
												Li.push_back(node);
											}
											contains = false;
										}
									}
								}
							}
						}
						F_Gf.clear();
					}
				}
			}
		}
	}

	ticksFinalFNM = clock() - ticksFNM;
	timeAllFNM = (double) ticksFinalFNM/CLOCKS_PER_SEC;
	//cout << "time in FindNewMapping: " << timeAllFNM << endl;

}




void UsePreviousMapping(unsigned int i, SOLUTION &solutioni, SOLUTION &solutioni_1, vector<TreeNode *> &Li, vector<TreeNode *> &Li_1, int window_size_increment)
{
	double timeAllUPM = 0;
	unsigned long ticksUPM = 0;
	unsigned long ticksFinalUPM = 0;

	ticksUPM = clock();

	SOL solutionif;
	MU mue_sol;

	vector<data_t> nodesHi;
	for(node_id nHi = 0; nHi < H[i].NodeCount(); nHi++)
	{
		nodesHi.push_back(H[i].GetNodeAttr(nHi));
	}


	for(unsigned int j = i - 1; j <= G.size() - H.size() + i - 1 + window_size_increment && j < G.size(); j++)  // search for previously found solutions within the window corresponding to H_{i-1}
	{	
		for(unsigned int k = 0; k < solutioni_1.size(); k++)
		{

			if(j == solutioni_1[k].first)
			{

				SOL solu = solutioni_1[k].second;                      // contains all previous mappings

				for(unsigned int a = 0; a < solu.size(); a++)
				{

					// iterate through all previous mappings

					MU mue = solu[a];

					// take nodes present in H[i]; filter using node attributes
					MU mue_new;
					for(unsigned int it = 0; it < mue.size(); it++)
					{
						auto iter = find(nodesHi.begin(), nodesHi.end(), mue[it].first);
						if(iter != nodesHi.end()){
							mue_new.push_back(make_pair(*iter, mue[it].second));
						}
					}


					//cout << "mue_new:" << endl;  		// mapping corresponding to G[j]; 
										// it's validity for G[f] is verified by finding the nodes 
										// present in mue_new as well as all edges between them in G[f] 

					vector<data_t> mue_new_first, mue_new_second;
					for(unsigned int i = 0; i < mue_new.size(); i++)
					{	
						mue_new_first.push_back(mue_new[i].first);
						mue_new_second.push_back(mue_new[i].second);
					}

					// check 
					// the 
					// future
					for(unsigned int f = j+1; f <= G.size() - H.size() + i + window_size_increment && f < G.size(); f++)     // sliding window to look into the future
					{
						// nodes of graph G[f]
						vector<data_t> nodesGf;
						for(node_id nGf = 0; nGf < G[f].NodeCount(); nGf++)
						{
							nodesGf.push_back(G[f].GetNodeAttr(nGf));
						}


						// look for nodes present in mue_new in graph G[f]
						vector<data_t> nodesfound_second;
						for(unsigned int l = 0; l < mue_new.size(); l++) 
						{
							auto it = find(nodesGf.begin(), nodesGf.end(), mue_new[l].second);   // here the search for a node in Gf is not in O(1) but in O(n)
							if(it != nodesGf.end()) 
							{
								nodesfound_second.push_back(mue_new[l].second);
							}
						}

						// look for edges between nodes present in mue_new in graph G[f]
						vector<vector<data_t> > neighsHi, neighsGf;
						vector<data_t> neighHi, neighGf;
						int mue_new_edgeCount_Hi = 0;
						int mue_new_edgeCount_Gf = 0;
						for(unsigned int k = 0; k < mue_new.size(); k++)
						{
							for(unsigned int l = 0; l < mue_new.size(); l++)
							{
								if(H[i].HasEdge(H[i].GetNodeId(mue_new[k].first), H[i].GetNodeId(mue_new[l].first)))
								{
									neighHi.push_back(mue_new[l].first);
								}
								if(G[f].HasEdge(G[f].GetNodeId(mue_new[k].second), G[f].GetNodeId(mue_new[l].second)))
								{
									neighGf.push_back(mue_new[l].second);
								}

							}

							mue_new_edgeCount_Hi += neighHi.size();    // number of edges between nodes present in mue_new in H[i] 
							mue_new_edgeCount_Gf += neighGf.size();    // number of edges between nodes present in mue_new in G[f] 

							if(!neighHi.empty())
								neighsHi.push_back(neighHi);
							if(!neighGf.empty())
								neighsGf.push_back(neighGf);

							neighHi.clear();
							neighGf.clear();
						}

						
						int index_neigh_Hi, index_neigh_Gf;
						vector<pair<data_t, data_t> > HiEdgesFound;

						for(unsigned int k = 0; k < neighsHi.size() && k < neighsGf.size(); k++)
						{
							for(unsigned int l = 0; l < neighsHi[k].size() && l < neighsGf[k].size(); l++)
							{
								auto iter_Hi = find(mue_new_first.begin(), mue_new_first.end(), neighsHi[k][l]);
								if(iter_Hi != mue_new_first.end())
								{
									index_neigh_Hi = distance(mue_new_first.begin(), iter_Hi);	
								}

								auto iter_Gf = find(mue_new_second.begin(), mue_new_second.end(), neighsGf[k][l]);
								if(iter_Gf != mue_new_second.end())
								{
									index_neigh_Gf = distance(mue_new_second.begin(), iter_Gf);	
								}

								if(index_neigh_Hi == index_neigh_Gf)
								{
									HiEdgesFound.push_back(make_pair(mue_new_first[k], neighsHi[k][l]));	
								}

							}
						}


						unsigned int it, c;
						if(nodesfound_second.size() == mue_new.size() && HiEdgesFound.size() == mue_new_edgeCount_Hi)
						{
							for(it=c= 0; it < nodesfound_second.size(), c < mue_new.size(); it++, c++)
							{
								node_id n = G[f].GetNodeId(nodesfound_second[it]);
								mue_sol.push_back(make_pair(mue_new[c].first, nodesfound_second[it]));
							}

							solutionif.push_back(mue_sol);
							mue_sol.clear();
							solutioni.push_back(make_pair(f,solutionif));


							bool contains = false;
							for(SOL::iterator itsol=solutionif.begin(); itsol!=solutionif.end();++itsol)
							{
								
								TreeNode *node = new TreeNode;
								
								for(auto it = Li_1.begin(); it != Li_1.end(); it++)
								{
									if((*it)->value == mue && (*it)->dg_index == j)               // to identify the correct parent node
									{                     
										node->parent = *it;                                   // use previous mapping
										node->value = *itsol;
										node->pg_index = i;
										node->dg_index = f;
										
										for(auto iter = (*it)->children.begin(); iter != (*it)->children.end(); ++iter)
										{
											if((*iter)->value == *itsol && (*iter)->dg_index == f)
											{
												contains = true;
											}

										}

										if(!contains)
										{
											(*it)->children.push_back(node);
											Li.push_back(node);
										}
										contains = false;
									}
								}
							}

							solutionif.clear();
						}
					}
				}
			}
		}
	}

	ticksFinalUPM = clock() - ticksUPM;
	timeAllUPM = (double) ticksFinalUPM/CLOCKS_PER_SEC;
	//cout << "time in UsePreviousMapping: " << timeAllUPM << endl;
}


void MapNextPatternGraphs(unsigned int i, SOLUTION &solutioni, SOLUTION &solutioni_1, vector<TreeNode *> &Li, vector<TreeNode *> &Li_1, int window_size_increment)
{

	//find common nodes between Hi and Hi-1 based on node attributes
	node_id nHi;
	node_id nHi_1;
	vector<node_id> NAMi;
	vector<node_id> NNMi;
	vector<node_id> nodesHi;

	for(nHi = 0; nHi < H[i].NodeCount(); nHi++) 
	{
		nodesHi.push_back(nHi);
		for(nHi_1 = 0; nHi_1 < H[i-1].NodeCount(); nHi_1++) 
		{
			if(H[i].GetNodeAttr(nHi) == H[i-1].GetNodeAttr(nHi_1)) 
			{
				NAMi.push_back(nHi);
			}
		}
	}

	for(auto itHi = nodesHi.begin(); itHi != nodesHi.end(); itHi++) 
	{
		if(find(NAMi.begin(), NAMi.end(), *itHi) == NAMi.end())               // if an element of Hi is not in NAMi, put it in NNMi
			NNMi.push_back(*itHi);
	}

	if(!NNMi.empty())
		FindNewMapping(i, NAMi, NNMi, solutioni, solutioni_1, Li, Li_1, window_size_increment);
	else
		UsePreviousMapping(i, solutioni, solutioni_1, Li, Li_1, window_size_increment);
}


void GetPathsRecur(vector<tuple<MU, unsigned int, unsigned int> > path, TreeNode* node, unsigned int pathLen, int &num){

	path.insert(path.end(), make_tuple(node->value, node->pg_index, node->dg_index));
	pathLen++;

	if(node->children.empty())
	{
		if(pathLen == H.size() + 1)
		{
			num++;
			for(unsigned int x = 1; x < path.size(); x++) 
			{
					for(unsigned int y = 0; y < (get<0>(path[x])).size(); y++)
					{
						cout << "H["<<get<1>(path[x])<<"]: " << (get<0>(path[x]))[y].first << " :: G["<<get<2>(path[x])<<"]: " << (get<0>(path[x]))[y].second << ", ";
					}
				cout << endl;
			}
		  	cout << endl;
		}
	}
	else 
	{
		for(unsigned int x = 0; x < node->children.size(); x++)
		{
			GetPathsRecur(path, node->children[x], pathLen, num);
		}
	}

}

int GetPaths(TreeNode *node){
	int numPaths = 0;
	vector<tuple<MU, unsigned int, unsigned int> > path;
	
	GetPathsRecur(path, node, 0, numPaths);

	return numPaths;
}


int main(int argc, char** argv){

	if(argc < 6)
	{
		cout << "Usage: [target graph parameters] [pattern graph parameters] [window_size_increment]" << endl;
		cout << "Usage: [window_size_increment] can be zero in which case it would be calculated using nb snapshots of target and pattern graphs" << endl;
		return -1;
	}

	string targetGraphFolderName, patternGraphFolderName;

	targetGraphFolderName = argv[1];
	int nbTimeStepsTarget = atoi(argv[2]);

	patternGraphFolderName = argv[3];
	int nbTimeStepsPattern = atoi(argv[4]);

	// for ablation study
	int window_size_increment = atoi(argv[5]);

	for(int a = 1; a <= nbTimeStepsTarget; a++)
	{	
			
		ifstream dg(targetGraphFolderName + "/dg" + to_string(a) + ".txt");
			
		StreamARGLoader<data_t, Empty> dg_loader(dg);                         // node attributes are strings (data_t) and edge attributes are empty (Empty)
		ARGraph<data_t, Empty> dg_graph(&dg_loader);

		G.push_back(dg_graph);
	}



	for(int a = 1; a <= nbTimeStepsPattern; a++)
	{
			
		ifstream pg(patternGraphFolderName + "/pg" + to_string(a) + ".txt");

		StreamARGLoader<data_t, Empty> pg_loader(pg);
		ARGraph<data_t, Empty> pg_graph(&pg_loader);

		H.push_back(pg_graph);
	}


	MU emptymu = {{"", ""}};

	TreeNode *root = new TreeNode;
	root->parent = NULL;
	root->value = emptymu;
	root->dg_index = -1;   // invalid index
	root->pg_index = -1;

	vector<vector<TreeNode *>> L = vector<vector<TreeNode *>> (H.size());
	vector<SOLUTION> S = vector<SOLUTION> (H.size());

	double timeAll = 0;
	unsigned long ticks = 0;
	unsigned long ticksFinal = 0;
	int nbSolutions;

	ticks = clock();                                  // total time except graph file reading

	for (unsigned int i = 0; i < H.size(); i++)
	{
		if(i == 0)
		{
			MapFirstPatternGraph(i, S[i], root, window_size_increment);
			L[i] = root->children;
		}

		if (i > 0) 
		{
			MapNextPatternGraphs(i, S[i], S[i-1], L[i], L[i-1], window_size_increment);
		}
	}

	nbSolutions = GetPaths(root);

	ticksFinal = clock() - ticks;
	timeAll =  (double) ticksFinal/CLOCKS_PER_SEC;

	cout << "no of solutions: " << nbSolutions << endl;
	cout << "time: " << timeAll << endl;

	//cout << timeAll << "," << nbSolutions << endl;

	G.clear();
	H.clear();

}
