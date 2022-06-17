//
//  vf3_sub_state.hpp
//  VF2Plus
//
//  Created by Vincenzo Carletti on 12/11/14.
//  Copyright (c) 2014 Vincenzo Carletti. All rights reserved.
//

#ifndef VF3_SUB_STATE_HPP
#define VF3_SUB_STATE_HPP

#include <cstring>
#include <iostream>
#include <vector>
#include "argraph.hpp"

typedef unsigned char node_dir_t;
#define NODE_DIR_NONE 0
#define NODE_DIR_IN	1
#define NODE_DIR_OUT 2
#define NODE_DIR_BOTH 3

using namespace std;

/*static void print_core(node_id* n1, node_id* n2, int size){
  for(int k = 0; k < size; k++){
    if(n1[k] != NULL_NODE )
      std::cout<<n2[n1[k]]<<","<<n1[k]<<":";
  }
}*/

/*----------------------------------------------------------
 * class VF3SubState
 * A representation of the SSR current state
 * See vf2_state.cc for more details.
 ---------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor = EqualityComparator<char, char>,
typename EdgeComparisonFunctor = EqualityComparator<Edge1, Edge2> >
//typename EdgeComparisonFunctor = EqualityComparator<string, string> >
class VF3SubState
{
private:
  //Comparison functors for nodes and edges
  NodeComparisonFunctor nf;
  EdgeComparisonFunctor ef;
  
  //Graphs to analyze
  ARGraph<Node1, Edge1> *g1;
  ARGraph<Node2, Edge2> *g2;
  
  const VF3SubState* parent;
  bool used;
  
  //Size of each graph
  int n1, n2;
  
  node_id *order;     //Order to traverse node on the first graph
  
  //CORE SET SIZES
  int core_len;       //Current length of the core set
  int orig_core_len;  //Core set length of the previous state
  int *core_len_c;    //Core set length for each class
  
  int added_node1;    //Last added node
  
  node_dir_t* dir;        //Node coming set. Direction into the terminal set.
  node_id* predecessors;  //Previous node in the ordered sequence connected to a node
  
  //TERMINAL SET SIZE
  //BE AWARE: Core nodes are also counted by these
  //GLOBAL SIZE
  int t2in_len, t2both_len, t2out_len; //Len of Terminal set for the second graph
  int *t1in_len, *t1both_len, *t1out_len; //Len of Terminal set for the first graph for each level
                                          //SIZE FOR EACH CLASS
  int *t2both_len_c, *t2in_len_c, *t2out_len_c;     //Len of Terminal set for the second graph for each class
  int **t1both_len_c, **t1in_len_c, **t1out_len_c;  //Len of Terminal set for the first graph for each class and level
  
  //Used for terminal set size evaluation
  int *termout2_c, *termin2_c, *new2_c;
  int **termout1_c, **termin1_c, **new1_c;
  int *termin1, *termout1, *new1;
  
  node_id *core_1;
  node_id *core_2;
  
  //Terminal sets of the second graph
  //TERM IN
  node_id *in_2;
  //TERMI OUT
  node_id *out_2;
  
  //Vector of sets used for searching the successors
  //Each class has its set
  int last_candidate_index;
  
  /* Structures for classes */
  int *class_1;       //Classes for nodes of the first graph
  int *class_2;       //Classes for nodes of the second graph
  int classes_count;  //Number of classes
  
  long *share_count;  //Count the number of instances sharing the common sets
  
  vector<node_id> nodesinmuHi, nodesinmuGf, F_Hi, F_Gf;  // to define mu for the initial state

  //PRIVATE METHODS
  void BackTrack();
  void ComputeFirstGraphTraversing();
  void ComputeGraphTraversing(vector<node_id> nodesinmuHi, vector<node_id> nodesinmuGf, vector<node_id> F_Hi, vector<node_id> F_Gf);
  void UpdateTerminalSetSize(node_id node, node_id level, bool* in_1, bool* out_1, bool* inserted);
  void UpdateTerminalSetSize(node_id node, node_id level, bool* in_1, bool* out_1, bool* inserted, vector<node_id> nodesinmuHi);
  void print_terminal(int c);
  
public:
  static long long instance_count;
  VF3SubState(ARGraph<Node1, Edge1> *g1, ARGraph<Node2, Edge2> *g2,
                int* class_1, int* class_2, int nclass,
                node_id* order = NULL);
  VF3SubState(ARGraph<Node1, Edge1> *g1, ARGraph<Node2, Edge2> *g2,
                  int* class_1, int* class_2, int nclass,
                  node_id* order, vector<node_id> nodesinmuHi, vector<node_id> nodesinmuGf, vector<node_id> F_Hi, vector<node_id> F_Gf);
  VF3SubState(const VF3SubState &state);
  ~VF3SubState();
  ARGraph<Node1, Edge1> *GetGraph1() { return g1; }
  ARGraph<Node2, Edge2> *GetGraph2() { return g2; }
  bool NextPair(node_id *pn1, node_id *pn2,
                node_id prev_n1=NULL_NODE, node_id prev_n2=NULL_NODE);
  bool IsFeasiblePair(node_id n1, node_id n2);
  void AddPair(node_id n1, node_id n2);
  bool IsGoal() { return core_len==n1; };
  bool IsDead();
  int CoreLen() { return core_len; }
  void GetCoreSet(node_id c1[], node_id c2[]);
  VF3SubState* GetParent();
  bool IsUsed(){return used;}
  void SetUsed(){used = true;}
};

/*----------------------------------------------------------
 * Methods of the class VF3SubState
 ---------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
void VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::print_terminal(int c){
  std::cout<<"\nClass: " << c << " Core_len: " << core_len_c[c];
  std::cout<<" t1both_len: " << t1both_len_c[core_len][c] << " t2both_len " << t2both_len_c[c];
  std::cout<<" t1out_len: " << t1out_len_c[core_len][c] << " t2out_len " << t2out_len_c[c];
  std::cout<<" t1in_len: " << t1in_len_c[core_len][c] << " t2in_len " << t2in_len_c[c];
  
}

/*----------------------------------------------------------
 * VF3SubState::VF3SubState(g1, g2)
 * Constructor. Makes an empty state.
 ---------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::VF3SubState(ARGraph<Node1, Edge1> *ag1, ARGraph<Node2, Edge2> *ag2, int* class_1, int* class_2, int nclass, node_id* order)
{ //cout << "in empty constructor" << endl;
  assert(class_1!=NULL && class_2!=NULL);
  
  VF3SubState::instance_count=1;
  g1=ag1;
  g2=ag2;
  n1=g1->NodeCount();
  n2=g2->NodeCount();
  last_candidate_index = 0;
  
  this->order = order;
  this->class_1 = class_1;
  this->class_2 = class_2;
  this->classes_count = nclass;
  parent = NULL;                                         // because empty state
  used = false;
  core_len=orig_core_len=0;                              //Current length of the core set = Core set length of the previous state
  t2both_len=t2in_len=t2out_len=0;
  
  //Creation of sets
  t1both_len = new int[n1+1];                             //Length of Terminal set for the first graph for each level
  t1in_len = new int[n1+1];
  t1out_len = new int[n1+1];
  
  termin1 = (int*) calloc(n1, sizeof(int));
  termout1 = (int*) calloc(n1, sizeof(int));
  new1 = (int*) calloc(n1, sizeof(int));
  
  t1both_len_c = (int**)malloc((n1+1)*sizeof(int*));             //Length of Terminal set for the first graph for each class end level
  t1in_len_c = (int**)malloc((n1+1)*sizeof(int*));
  t1out_len_c = (int**)malloc((n1+1)*sizeof(int*));
  
  termin1_c = (int**)malloc(n1*sizeof(int*));
  termout1_c = (int**)malloc(n1*sizeof(int*));
  new1_c = (int**)malloc(n1*sizeof(int*));
  
  core_len_c = (int*) calloc(classes_count, sizeof(int));                   //Core set length for each class
  t2both_len_c = (int*) calloc(classes_count, sizeof(int));
  t2in_len_c = (int*) calloc(classes_count, sizeof(int));
  t2out_len_c = (int*) calloc(classes_count, sizeof(int));
  termout2_c = new int[classes_count];
  termin2_c = new int[classes_count];
  new2_c = new int[classes_count];
  
  added_node1=NULL_NODE;
  
  core_1=new node_id[n1];                          //core set for first graph
  core_2=new node_id[n2];                          // core set for second graph
  in_2=new node_id[n2];                            // Terminal sets of the second graph: TERM IN
  out_2=new node_id[n2];                           // Terminal sets of the second graph: TERM OUT
  dir = new node_dir_t[n1];
  predecessors = new node_id[n1];
  share_count = new long;
  
  int i;
  for(i=0; i<=n1; i++)
    {
    if(i<n1){
      core_1[i]=NULL_NODE;
        termin1_c[i] = (int*) calloc(classes_count, sizeof(int));
        termout1_c[i] = (int*) calloc(classes_count, sizeof(int));
        new1_c[i] = (int*) calloc(classes_count, sizeof(int));
    }
      t1both_len_c[i] = (int*) calloc(classes_count, sizeof(int));
      t1in_len_c[i] = (int*) calloc(classes_count, sizeof(int));
      t1out_len_c[i] = (int*) calloc(classes_count, sizeof(int));
    }
  
  for(i=0; i<n2; i++)
    {
    core_2[i]=NULL_NODE;
    in_2[i]=0;
    out_2[i]=0;
    }
  
  ComputeFirstGraphTraversing();      // preprocess pattern graph
  *share_count = 1;
}

/*
 * VF3SubState constructor with nodesinmuHi and nodesinmuGf
 *
 */

template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::VF3SubState(ARGraph<Node1, Edge1> *ag1, ARGraph<Node2, Edge2> *ag2,
																				int* class_1, int* class_2, int nclass, node_id* order,
																				vector<node_id> nodesinmuHi, vector<node_id> nodesinmuGf,
																				vector<node_id> F_Hi, vector<node_id> F_Gf){
   //cout << "in new constructor" << endl;
   assert(class_1!=NULL && class_2!=NULL);

   VF3SubState::instance_count=1;                         // to make a new SSR
   g1=ag1;
   g2=ag2;
   n1=g1->NodeCount();
   n2=g2->NodeCount();

   last_candidate_index = 0;         // TODO

   this->order = order;
   this->class_1 = class_1;
   this->class_2 = class_2;
   this->classes_count = nclass;
   parent = NULL;                                         // because initial state
   used = false;                     // TODO
   orig_core_len=0;                                   //Core set length of the previous state

   t2both_len=t2in_len=t2out_len=0;   //TODO same as feasibility set?

   core_len=0;              // current length of the core set

   //Creation of sets
   t1both_len = new int[n1+1];                             //Length of Terminal set for the first graph for each level
   t1in_len = new int[n1+1];
   t1out_len = new int[n1+1];

   termin1 = (int*) calloc(n1, sizeof(int));
   termout1 = (int*) calloc(n1, sizeof(int));
   new1 = (int*) calloc(n1, sizeof(int));

   t1both_len_c = (int**)malloc((n1+1)*sizeof(int*));             //Length of Terminal set for the first graph for each class end level
   t1in_len_c = (int**)malloc((n1+1)*sizeof(int*));
   t1out_len_c = (int**)malloc((n1+1)*sizeof(int*));

   termin1_c = (int**)malloc(n1*sizeof(int*));
   termout1_c = (int**)malloc(n1*sizeof(int*));
   new1_c = (int**)malloc(n1*sizeof(int*));

   core_len_c = (int*) calloc(classes_count, sizeof(int));                   //Core set length for each class
   t2both_len_c = (int*) calloc(classes_count, sizeof(int));
   t2in_len_c = (int*) calloc(classes_count, sizeof(int));
   t2out_len_c = (int*) calloc(classes_count, sizeof(int));
   termout2_c = new int[classes_count];
   termin2_c = new int[classes_count];
   new2_c = new int[classes_count];

   added_node1=NULL_NODE;

   core_1=new node_id[n1];                          //core set for first graph
   core_2=new node_id[n2];                          // core set for second graph
   in_2=new node_id[n2];                            // Terminal sets of the second graph: TERM IN
   out_2=new node_id[n2];                           // Terminal sets of the second graph: TERM OUT
   dir = new node_dir_t[n1];
   predecessors = new node_id[n1];
   share_count = new long;

   int i;
   for(i=0; i<=n1; i++)
   {
	 if(i<n1){
		 core_1[i]=NULL_NODE;
		 termin1_c[i] = (int*) calloc(classes_count, sizeof(int));
         termout1_c[i] = (int*) calloc(classes_count, sizeof(int));
         new1_c[i] = (int*) calloc(classes_count, sizeof(int));
     }
     t1both_len_c[i] = (int*) calloc(classes_count, sizeof(int));
     t1in_len_c[i] = (int*) calloc(classes_count, sizeof(int));
     t1out_len_c[i] = (int*) calloc(classes_count, sizeof(int));
   }

   for(i=0; i<n2; i++)
   {
	 core_2[i]=NULL_NODE;
     in_2[i]=0;
     out_2[i]=0;
   }



   ComputeGraphTraversing(nodesinmuHi, nodesinmuGf, F_Hi, F_Gf);                // preprocess pattern graph and compute sets for target graph


   *share_count = 1;



}




/*----------------------------------------------------------
 * VF3SubState::VF3SubState(state)
 * Copy constructor.
 ---------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::VF3SubState(const VF3SubState &state)
{ //cout << "in copy constructor" << endl;
  g1=state.g1;
  g2=state.g2;
  n1=state.n1;
  n2=state.n2;
  
  order=state.order;
  class_1 = state.class_1;
  class_2 = state.class_2;
  classes_count = state.classes_count;
  parent = &state;                                               // previous state becomes the parent of the current state
  used = false;
  VF3SubState::instance_count++;
  
  last_candidate_index = state.last_candidate_index;
  
  core_len=orig_core_len=state.core_len;
  //cout << "core_len: " << core_len << endl;
  
  t1in_len=state.t1in_len;
  t1out_len=state.t1out_len;
  t1both_len=state.t1both_len;
  
  t2in_len=state.t2in_len;
  t2out_len=state.t2out_len;
  t2both_len=state.t2both_len;
  
  core_len_c = state.core_len_c;
  t1both_len_c = state.t1both_len_c;
  t2both_len_c = state.t2both_len_c;
  t1in_len_c = state.t1in_len_c;
  t2in_len_c = state.t2in_len_c;
  t1out_len_c = state.t1out_len_c;
  t2out_len_c = state.t2out_len_c;
  
  termout1_c = state.termout1_c;
  termout2_c = state.termout2_c;
  termin1_c = state.termin1_c;
  termin2_c = state.termin2_c;
  new1_c = state.new1_c;
  new2_c = state.new2_c;
  
  termin1 = state.termin1;
  termout1 = state.termout1;
  new1 = state.new1;
  
  added_node1=NULL_NODE;
  
  core_1=state.core_1;
  core_2=state.core_2;
  in_2=state.in_2;
  out_2=state.out_2;
  //cout << "core_1: " << *core_1<<", core_2: "<<*core_2<<", in_2: " <<*in_2<<", out_2: "<<*out_2 << endl;
  dir = state.dir;
  predecessors = state.predecessors;
  share_count=state.share_count;
  
  ++ *share_count;
  
  //cout << "share_count: " << *share_count << endl;
  //cout << "instance_count: " << instance_count << endl;

}


/*---------------------------------------------------------------
 * VF3SubState::~VF3SubState()
 * Destructor.
 --------------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::~VF3SubState()
{
	//cout << "in destructor" << endl;

	//cout << "share_count: " << *share_count << endl;
	 //cout << "instance_count: " << instance_count << endl;
  
  if(-- *share_count > 0){               // when destructing non empty state
	  //cout << "share_count reduced: " << *share_count << endl;
	  BackTrack();
  }
  
  if (*share_count == 0)                // when destructing empty state
    { delete [] core_1;
      delete [] core_2;
      delete [] in_2;
      delete [] out_2;
      delete [] dir;
      delete [] predecessors;
      delete [] t1both_len;
      delete [] t1in_len;
      delete [] t1out_len;
      delete [] termin1;
      delete [] termout1;
      delete [] new1;
      
      for(int i = 0; i <= n1; i++){
        delete [] t1both_len_c[i];
        delete [] t1in_len_c[i];
        delete [] t1out_len_c[i];
        if(i< n1){
          delete [] termin1_c[i];
          delete [] termout1_c[i];
          delete [] new1_c[i];
        }
      }
      
      delete [] t1both_len_c;
      delete [] t1in_len_c;
      delete [] t1out_len_c;
      delete [] termin1_c;
      delete [] termout1_c;
      delete [] new1_c;
      delete [] t2both_len_c;
      delete [] t2in_len_c;
      delete [] t2out_len_c;
      delete [] core_len_c;
      delete [] termin2_c;
      delete [] termout2_c;
      delete [] new2_c;
      
      delete share_count;
    }
}

template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
void VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,
EdgeComparisonFunctor>::UpdateTerminalSetSize(node_id node, node_id level, bool* in_1, bool* out_1, bool* inserted){
  node_id i, neigh, c_neigh;
  node_id in1_count, out1_count;
  
  //Updating Terminal set size count And degree
  in1_count = g1->InEdgeCount(node);
  out1_count = g1->OutEdgeCount(node);

  //Updating Inner Nodes not yet inserted
  for (i = 0; i < in1_count; i++)
    {
    //Getting Neighborhood
    neigh = g1->GetInEdge(node,i);
    c_neigh = class_1[neigh];
    
    if(!inserted[neigh])
      {
      if (in_1[neigh]){
        termin1[level]++;                           // just change the size
        termin1_c[level][c_neigh]++;
      }
      if (out_1[neigh]){
        termout1[level]++;
        termout1_c[level][c_neigh]++;
      }
      if (!in_1[neigh] && !out_1[neigh]){
        new1[level]++;
        new1_c[level][c_neigh]++;
      }
      }
    }
  
  //Updating Outer Nodes not yet inserted
  for (i = 0; i < out1_count; i++)
    {
    //Getting Neighborhood
    neigh = g1->GetOutEdge(node,i);
    c_neigh = class_1[neigh];
    if(!inserted[neigh])
      {
      if (in_1[neigh]){
        termin1[level]++;
        termin1_c[level][c_neigh]++;
      }
      if (out_1[neigh]){
        termout1[level]++;
        termout1_c[level][c_neigh]++;
      }
      if (!in_1[neigh] && !out_1[neigh]){
        new1[level]++;
        new1_c[level][c_neigh]++;
      }
      }
    }
}



//Try to have predetermined in1 and out1, without having to calculate it at each iteration
//Their size at each level of the search tree is predetermined
//In this way I need only know the order of choice and the size of in1 and out1
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
void VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::ComputeFirstGraphTraversing(){
  //The algorithm start with the node with the maximum degree
  node_id depth, i;
  node_id node;	//Current Node
  int node_c; //Class of the current node
  bool* inserted = new bool[n1];
  bool *in, *out; //Internal Terminal Set used for updating the size of
  in = new bool[n1];
  out = new bool[n1];
  
  //Init vectors and variables
  node = 0;
  node_c = 0;
  
  t1in_len[0] = 0;
  t1out_len[0] = 0;
  t1both_len[0] = 0;
  
  for(i = 0; i < n1; i++)
    {
    in[i] = false;
    out[i] = false;
    dir[i] = NODE_DIR_NONE;
    inserted[i] = false;
    predecessors[i] = NULL_NODE;
    }
  
  /* Following the imposed node order */
  for(depth = 0; depth < n1; depth++)
    {
	  //cout << "depth: " << depth << endl;
    node = order[depth];
    node_c = class_1[node];
    inserted[node] = true;                            // inserted where ??? in terminal set
    
    UpdateTerminalSetSize(node, depth, in, out, inserted);
    
    //Updating counters for next step
    t1in_len[depth +1] = t1in_len[depth];
    t1out_len[depth +1] = t1out_len[depth];
    t1both_len[depth +1] = t1both_len[depth];
    for (int j = 0; j < classes_count; j++)
      {
      t1in_len_c[depth +1][j] = t1in_len_c[depth][j];
      t1out_len_c[depth +1][j] = t1out_len_c[depth][j];
      t1both_len_c[depth +1][j] = t1both_len_c[depth][j];
      }
    //Inserting the node
    //Terminal set sizes depends on the depth
    // < depth non sono nell'insieme
    // >= depth sono nell'insieme
    if (!in[node])
      {
      in[node]=true;
      t1in_len[depth+1]++;
      t1in_len_c[depth+1][node_c]++;
      if (out[node]){
        t1both_len[depth+1]++;
        t1both_len_c[depth+1][node_c]++;
      }
      }
    
    if (!out[node])
      {
      out[node]=true;
      t1out_len[depth+1]++;
      t1out_len_c[depth+1][node_c]++;
      if (in[node]){
        t1both_len[depth+1]++;
        t1both_len_c[depth+1][node_c]++;
      }
      }
    
    //Updating terminal sets
    int i, other, other_c;
    other_c = -1;
    for(i=0; i<g1->InEdgeCount(node); i++)
      {
      other=g1->GetInEdge(node, i);
      if (!in[other])
        {
        other_c = class_1[other];
        in[other]=true;
        t1in_len[depth+1]++;
        t1in_len_c[depth+1][other_c]++;
        if(!inserted[other])
        {
          if(predecessors[other] == NULL_NODE)
          {
            dir[other] = NODE_DIR_IN;
            predecessors[other] = node;
          }
        }
        if (out[other]){
          t1both_len[depth+1]++;
          t1both_len_c[depth+1][other_c]++;
        }
        }
      }
    
    for(i=0; i<g1->OutEdgeCount(node); i++)
      {
      other=g1->GetOutEdge(node, i);
      if (!out[other])
        {
        other_c = class_1[other];
        out[other]=true;
        t1out_len[depth+1]++;
        t1out_len_c[depth+1][other_c]++;
        if(!inserted[other])
        {
          if(predecessors[other] == NULL_NODE)
          {
            predecessors[other] = node;
            dir[other] = NODE_DIR_OUT;
          }
        }
        if (in[other]){
          t1both_len[depth+1]++;
          t1both_len_c[depth+1][other_c]++;
        }
        }
      }

    /*cout << "core_len: " << core_len << endl;
    cout <<"t1both_len["<<core_len<<"]: " <<t1both_len[core_len] << ", t2both_len: " << t2both_len << endl;
    cout <<"t1out_len["<<core_len<<"]: " <<t1out_len[core_len] << ", t2out_len: " << t2out_len << endl;
    cout <<"t1in_len["<<core_len<<"]: " <<t1in_len[core_len] << ", t2in_len: " << t2in_len << endl;*/
    }
  
  delete [] in;
  delete [] out;
  delete [] inserted;
}


template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
void VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,
EdgeComparisonFunctor>::UpdateTerminalSetSize(node_id node, node_id level, bool* in_1, bool* out_1, bool* inserted,
												vector<node_id> nodesinmuHi){
  node_id i, neigh, c_neigh;
  node_id in1_count, out1_count;

  //Updating Terminal set size count And degree
  in1_count = g1->InEdgeCount(node);
  out1_count = g1->OutEdgeCount(node);

  //cout <<"UpdateTerminalSetSize: " << "in1_count: " << in1_count << ", out1_count: " << out1_count << " of node: " << node << endl;

  //Updating Inner Nodes not yet inserted
  for (i = 0; i < in1_count; i++)
    {
    //Getting Neighborhood
    neigh = g1->GetInEdge(node,i);
    c_neigh = class_1[neigh];

    //cout << "inner nodes" << endl;
    //cout << "neigh: " << neigh << endl;
    //cout<<"in_1["<<neigh<<"]: " <<in_1[neigh] << endl;
    //cout<<"out_1["<<neigh<<"]: " <<out_1[neigh] << endl;

    //auto it = find(nodesinmuHi.begin(), nodesinmuHi.end(), neigh);

    if(!inserted[neigh])// && it==nodesinmuHi.end())               // neigh not in core set of first graph
      {
    	//cout << "here_inner" << endl;
      if (in_1[neigh]){
        termin1[level]++;
        termin1_c[level][c_neigh]++;
      }
      if (out_1[neigh]){
        termout1[level]++;
        termout1_c[level][c_neigh]++;
      }
      if (!in_1[neigh] && !out_1[neigh]){
        new1[level]++;
        new1_c[level][c_neigh]++;

      }
      }
    }

  /*cout << "termin1["<<level<<"]: " << termin1[level] << endl;
  cout << "termout1["<<level<<"]: " << termout1[level] << endl;
  cout << "new1["<<level<<"]: " << new1[level] << endl;*/

  //Updating Outer Nodes not yet inserted
  for (i = 0; i < out1_count; i++)
    {
    //Getting Neighborhood
    neigh = g1->GetOutEdge(node,i);
    c_neigh = class_1[neigh];

    //cout << "outer nodes" << endl;
    //cout << "neigh: " << neigh << endl;
    //cout<<"in_1["<<neigh<<"]: " <<in_1[neigh] << endl;
    //cout<<"out_1["<<neigh<<"]: " <<out_1[neigh] << endl;

    //auto it = find(nodesinmuHi.begin(), nodesinmuHi.end(), neigh);

    if(!inserted[neigh])// && it==nodesinmuHi.end())                 // neigh not in core set of first graph
      { //cout << "here_outer" << endl;
      if (in_1[neigh]){
        termin1[level]++;
        termin1_c[level][c_neigh]++;
      }
      if (out_1[neigh]){
        termout1[level]++;
        termout1_c[level][c_neigh]++;
      }
      if (!in_1[neigh] && !out_1[neigh]){
        new1[level]++;
        new1_c[level][c_neigh]++;
      }
      }
    }

  /*cout << "termin1["<<level<<"]: " << termin1[level] << endl;
  cout << "termout1["<<level<<"]: " << termout1[level] << endl;
  cout << "new1["<<level<<"]: " << new1[level] << endl;*/

}



/*
 * ComputeGraphTraversing: preprocess pattern graph: compute all sets for pattern graph
 * Also compute all sets for target graph until depth > nodesinmuHi.size()
 *
 */
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
void VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::ComputeGraphTraversing(vector<node_id> nodesinmuHi,
															vector<node_id> nodesinmuGf, vector<node_id> F_Hi, vector<node_id> F_Gf){

	node_id depth, i;
	node_id node1;	//Current Node
	node_id node2;  //node paired to node1
	int node1_c; //Class of the current node
	bool* inserted = new bool[n1];         // for node inserted in core set of first graph
	bool *in, *out;              //Internal Terminal Set used for updating the size of
	in = new bool[n1];
	out = new bool[n1];

	//Init vectors and variables
	node1 = 0;
	node1_c = 0;


	for(i = 0; i < n1; i++)
	{
	  in[i] = false;
	  out[i] = false;
	  dir[i] = NODE_DIR_NONE;
	  inserted[i] = false;
	  predecessors[i] = NULL_NODE;
	}

	/*
	 * for depth = 1
	 * take core set of first graph from nodesinmuHi
	 * core set of second graph from nodesinmuGf
	 *
	 * take terminal set of first graph from F_Hi
	 * terminal set of second graph from F_Gf
	 *
	 */

	core_len = nodesinmuHi.size();
	unsigned int a,b;
	for(a=b= 0; a < nodesinmuHi.size(), b < nodesinmuGf.size(); a++, b++)
	{

		node1 = nodesinmuHi[a];
		node2 = nodesinmuGf[b];
		node1_c = class_1[node1];
		core_len_c[node1_c]++;
		inserted[node1] = true;

		//Inserting nodes into the core set
		core_1[node1]=node2;
		core_2[node2]=node1;

		//cout << "core_1["<<node1<<"] = " <<core_1[node1] << endl;
		//cout << "core_2["<<node2<<"] = " <<core_2[node2] << endl;
    }



	for(unsigned int i = 0; i < nodesinmuHi.size(); i++){
		for(int j = 0; j < g1->OutEdgeCount(nodesinmuHi[i]); j++){
			node_id n = g1->GetOutEdge(nodesinmuHi[i], j);
			if(!inserted[n]){ // not in core set
				predecessors[nodesinmuHi[i]] = n;               // predecessor defined for nodes in core set
				//predecessors[n] = nodesinmuHi[i];
			}
		}
	}

	//t1in, t1out; in, out
	t1in_len[1] = F_Hi.size();
	t1out_len[1] = F_Hi.size();
	t1both_len[1] = F_Hi.size();
	for(unsigned int i = 0; i < F_Hi.size(); i++){

		int n1_c = class_1[F_Hi[i]];

		//cout << "n1_c = " << n1_c << endl;

		in[F_Hi[i]] = true;
		t1in_len_c[1][n1_c]++;
		out[F_Hi[i]] = true;
		t1out_len_c[1][n1_c]++;
		t1both_len_c[1][n1_c]++;
	}


	//cout << "core_1[0] === " << core_1[0] << endl;

	//t2in, t2out: in_2, out_2
	t2in_len = F_Gf.size();
	t2out_len = F_Gf.size();
	t2both_len = F_Gf.size();
	for(unsigned int j = 0; j < F_Gf.size(); j++){

		int n2_c = class_2[F_Gf[j]];

		//cout << "n2_c = " << n2_c << endl;

		in_2[F_Gf[j]] = true;
		t2in_len_c[n2_c]++;
		out_2[F_Gf[j]] = true;
		t2out_len_c[n2_c]++;
		t2both_len_c[n2_c]++;
	}

	//cout << "core_1[0] === " << core_1[0] << endl;

	/*cout << "depth: 0" << endl;
	cout <<"t1both_len["<<0<<"]: " <<t1both_len[0] << ", t2both_len: " << t2both_len << endl;
	cout <<"t1out_len["<<0<<"]: " <<t1out_len[0] << ", t2out_len: " << t2out_len << endl;
	cout <<"t1in_len["<<0<<"]: " <<t1in_len[0] << ", t2in_len: " << t2in_len << endl;*/


	/*
	 * for depth > 0
	 * take core set of first graph from NG1
	 * core set of second graph computed by VF3
	 *
	 * terminal set of first graph precomputed before matching  (like in ComputeFirstGraphTraversing)
	 * terminal set of second graph computed by VF3
	 *
	 */

	//t1in_len[1] = 0;
	//t1out_len[1] = 0;                // for depth = 0 is given by nodesinmuHi.size()
	//t1both_len[1] = 0;

	//cout << "core_1[0] === " << core_1[0] << endl;

	  /* Following the imposed node order */
	  //for(depth = 1; depth <= n1-core_len; depth++)
	    for(depth = core_len; depth < n1; depth++)
	    {
		  //cout << "depth: " << depth << endl;
	      node1 = order[depth];                       // order padded with zeros upto index < core_len
	      //cout << "node1: " << node1 << endl;
	      //node1 = order[depth-1];   // order starts from index 0
	      node1_c = class_1[node1];
	      inserted[node1] = true;

	    UpdateTerminalSetSize(node1, depth, in, out, inserted, nodesinmuHi);     // terminal set for current level

	    /*cout <<"t1both_len["<<depth<<"]: " <<t1both_len[depth] << ", t2both_len: " << t2both_len << endl;
	    cout <<"t1out_len["<<depth<<"]: " <<t1out_len[depth] << ", t2out_len: " << t2out_len << endl;
	    cout <<"t1in_len["<<depth<<"]: " <<t1in_len[depth] << ", t2in_len: " << t2in_len << endl;*/

	    if(depth == core_len){
		    t1in_len[depth] = t1in_len[1];
		    t1out_len[depth] = t1out_len[1];
		    t1both_len[depth] = t1both_len[1];
		    for (int j = 0; j < classes_count; j++)
		    {
		    	t1in_len_c[depth +1][j] = t1in_len_c[1][j];
		    	t1out_len_c[depth +1][j] = t1out_len_c[1][j];
		    	t1both_len_c[depth +1][j] = t1both_len_c[1][j];

		    	/*cout <<"t1both_len_c["<<depth+1<<"]["<<j<<"]: " <<t1both_len_c[depth +1][j] << endl;
		    	  cout <<"t1out_len_c["<<depth+1<<"]["<<j<<"]: " <<t1out_len_c[depth +1][j] << endl;
		    	  cout <<"t1in_len_c["<<depth+1<<"]["<<j<<"]: " <<t1in_len_c[depth +1][j] << endl;*/

		    }
	    }


	    //Updating counters for next step
	    t1in_len[depth +1] = t1in_len[depth];
	    t1out_len[depth +1] = t1out_len[depth];
	    t1both_len[depth +1] = t1both_len[depth];
	    for (int j = 0; j < classes_count; j++)
	      {
	      t1in_len_c[depth +1][j] = t1in_len_c[depth][j];
	      t1out_len_c[depth +1][j] = t1out_len_c[depth][j];
	      t1both_len_c[depth +1][j] = t1both_len_c[depth][j];

	      /*cout <<"t1both_len_c["<<depth+1<<"]["<<j<<"]: " <<t1both_len_c[depth +1][j] << endl;
	      cout <<"t1out_len_c["<<depth+1<<"]["<<j<<"]: " <<t1out_len_c[depth +1][j] << endl;
	      cout <<"t1in_len_c["<<depth+1<<"]["<<j<<"]: " <<t1in_len_c[depth +1][j] << endl;*/

	      }


	    //Inserting the node
	    //Terminal set sizes depends on the depth
	    // < depth non sono nell'insieme
	    // >= depth sono nell'insieme
	    if (!in[node1])
	      {
	      in[node1]=true;
	      t1in_len[depth+1]++;
	      t1in_len_c[depth+1][node1_c]++;
	      if (out[node1]){
	        t1both_len[depth+1]++;
	        t1both_len_c[depth+1][node1_c]++;
	      }
	      }

	    if (!out[node1])
	      {
	      out[node1]=true;
	      t1out_len[depth+1]++;
	      t1out_len_c[depth+1][node1_c]++;
	      if (in[node1]){
	        t1both_len[depth+1]++;
	        t1both_len_c[depth+1][node1_c]++;
	      }
	      }

	    //Updating terminal sets
	    int i, other, other_c;
	    other_c = -1;
	    for(i=0; i<g1->InEdgeCount(node1); i++)
	      {
	      other=g1->GetInEdge(node1, i);
	      //cout << "other: " << other << endl;
	      //cout <<"predecessors["<<other<<"]: " <<predecessors[other] << endl;

	      if (!in[other] && !inserted[other])         // if neighbor in core set don't put it in terminal set
	      //if (!in[other])
	        { //cout << "!in[other]" << endl;
	        other_c = class_1[other];
	        in[other]=true;
	        t1in_len[depth+1]++;
	        t1in_len_c[depth+1][other_c]++;
	        if(!inserted[other])
	        {
	          if(predecessors[other] == NULL_NODE)
	          {
	            dir[other] = NODE_DIR_IN;
	            predecessors[other] = node1;       // shouldn't be predecessors[node1] = other ???
	          }
	        }
	        if (out[other]){
	          t1both_len[depth+1]++;
	          t1both_len_c[depth+1][other_c]++;
	        }
	        }
	      }


	    /*cout << "first one" << endl;
	    for (int j = 0; j < classes_count; j++)
	    	    {
	    	    	cout <<"t1both_len_c["<<depth+1<<"]["<<j<<"]: " <<t1both_len_c[depth +1][j] << ", t2both_len_c: " << t2both_len_c[j] << endl;
	    	    	cout <<"t1out_len_c["<<depth+1<<"]["<<j<<"]: " <<t1out_len_c[depth +1][j] << ", t2out_len_c: " << t2out_len_c[j] << endl;
	    	    	cout <<"t1in_len_c["<<depth+1<<"]["<<j<<"]: " <<t1in_len_c[depth +1][j] << ", t2in_len_c: " << t2in_len_c[j] << endl;
	    	    }*/

	    for(i=0; i<g1->OutEdgeCount(node1); i++)
	      {
	      other=g1->GetOutEdge(node1, i);
	      if (!out[other] && !inserted[other])           // if neighbor in core set don't put it in terminal set
	      //if (!out[other])
	        {//cout << "!out[other]" << endl;
	        other_c = class_1[other];
	        out[other]=true;
	        t1out_len[depth+1]++;
	        t1out_len_c[depth+1][other_c]++;
	        if(!inserted[other])
	        {
	          if(predecessors[other] == NULL_NODE)
	          {
	            predecessors[other] = node1;
	            dir[other] = NODE_DIR_OUT;
	          }
	        }
	        if (in[other]){
	          t1both_len[depth+1]++;
	          t1both_len_c[depth+1][other_c]++;
	        }
	        }
	      }
	    /*cout << "next one" << endl;
	    for (int j = 0; j < classes_count; j++)
	    {
	    	cout <<"t1both_len_c["<<depth+1<<"]["<<j<<"]: " <<t1both_len_c[depth +1][j] << ", t2both_len_c: " << t2both_len_c[j] << endl;
	    	cout <<"t1out_len_c["<<depth+1<<"]["<<j<<"]: " <<t1out_len_c[depth +1][j] << ", t2out_len_c: " << t2out_len_c[j] << endl;
	    	cout <<"t1in_len_c["<<depth+1<<"]["<<j<<"]: " <<t1in_len_c[depth +1][j] << ", t2in_len_c: " << t2in_len_c[j] << endl;
	    }*/

	    /*cout << "depth: " << depth << endl;
	    cout << "core_len: " << core_len << endl;
	    cout <<"t1both_len["<<depth<<"]: " <<t1both_len[depth] << ", t2both_len: " << t2both_len << endl;
	   	cout <<"t1out_len["<<depth<<"]: " <<t1out_len[depth] << ", t2out_len: " << t2out_len << endl;
	   	cout <<"t1in_len["<<depth<<"]: " <<t1in_len[depth] << ", t2in_len: " << t2in_len << endl;*/
	 }
	  //cout << "core_1[0] === " << core_1[0] << endl;

	  /*cout << "depth: " << depth << endl;
	  cout << "core_len: " << core_len << endl;
	  cout <<"t1both_len["<<core_len<<"]: " <<t1both_len[core_len] << ", t2both_len: " << t2both_len << endl;
	  cout <<"t1out_len["<<core_len<<"]: " <<t1out_len[core_len] << ", t2out_len: " << t2out_len << endl;
	  cout <<"t1in_len["<<core_len<<"]: " <<t1in_len[core_len] << ", t2in_len: " << t2in_len << endl;*/


	  delete [] in;
	  delete [] out;
	  delete [] inserted;
}





template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
bool VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::NextPair(node_id *pn1, node_id *pn2,node_id prev_n1, node_id prev_n2)
{
  //cout << "in next pair" << endl;
  node_id curr_n1;
  node_id pred_pair; //Node mapped with the predecessor
  node_id pred_set_size = 0;
  int c = 0;
  pred_pair = NULL_NODE;
  

  //core_len indicates the depth of the research
  //cout << "core_len: " << core_len << endl;
  curr_n1 = order[core_len];
  c = class_1[curr_n1];
  
  //cout << "curr_n1: " << curr_n1 << endl;
  //cout <<"predecessors["<<curr_n1<<"]: " <<predecessors[curr_n1] << endl;
  //cout << "prev_n1: " << prev_n1 << ", prev_n2: " << prev_n2 << endl;

  //cout << "core_1[0] === " << core_1[0] << endl;

  if(predecessors[curr_n1] != NULL_NODE)
  {
	  if (prev_n2 == NULL_NODE)
		  last_candidate_index = 0;
	  else
	  {
		  last_candidate_index++; //Next Element
	  }

	  //cout << "last_candidate_index: " << last_candidate_index << endl;

	  pred_pair = core_1[predecessors[curr_n1]];
	  //cout << "pred_pair: " << pred_pair << endl;

	  //cout << "dir[curr_n1]: " << dir[curr_n1] << endl;

	  switch (dir[curr_n1])
      {
        	case NODE_DIR_IN:
        		pred_set_size = g2->InEdgeCount(pred_pair);
        		while(last_candidate_index < pred_set_size)
        		{
        			prev_n2 = g2->GetInEdge(pred_pair,last_candidate_index);
        			if(core_2[prev_n2] != NULL_NODE || class_2[prev_n2] != c)
        				last_candidate_index++;
        			else
        				break;
        		}
        	break;
        
        	case NODE_DIR_OUT:
        		pred_set_size = g2->OutEdgeCount(pred_pair);
        		while(last_candidate_index < pred_set_size)
        		{
        			prev_n2 = g2->GetOutEdge(pred_pair,last_candidate_index);
        			if(core_2[prev_n2] != NULL_NODE || class_2[prev_n2] != c)
        				last_candidate_index++;
        			else
        				break;
        		}
        	break;
      }
	  if(last_candidate_index >= pred_set_size)
		  return false;
  }
  else   // predecessors[curr_n1] == null
  {
	  //Recupero il nodo dell'esterno
	  if(prev_n2 == NULL_NODE)
		  prev_n2 = 0;
	  else
		  prev_n2++;
    
	  while (prev_n2<n2 && (core_2[prev_n2]!=NULL_NODE || class_2[prev_n2] != c))
      {
		  prev_n2++;
      }
  }
  //std::cout<<curr_n1 << " " << prev_n2 << " \n";
  
  //cout << "core_1[0] === " << core_1[0] << endl;

  if (prev_n2 < n2) {
    *pn1 = curr_n1;
    *pn2 = prev_n2;
    //cout<<"NP END: " <<curr_n1<<" " << prev_n2 << "\n" ;
    return true;
  }
  
  return false;
}


/*---------------------------------------------------------------
 * bool VF3SubState::IsFeasiblePair(node1, node2)
 * Returns true if (node1, node2) can be added to the state
 * NOTE:
 *   The attribute compatibility check (methods CompatibleNode
 *   and CompatibleEdge of ARGraph) is always performed
 *   applying the method to g1, and passing the attribute of
 *   g1 as first argument, and the attribute of g2 as second
 *   argument. This may be important if the compatibility
 *   criterion is not symmetric.
 --------------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
bool VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::IsFeasiblePair(node_id node1, node_id node2)
{
  //std::cout<<"\nIF: " <<node1<<" " << node2;
  //print_core(core_1, core_2, core_len);
  assert(node1<n1);
  assert(node2<n2);
  assert(core_1[node1]==NULL_NODE);
  assert(core_2[node2]==NULL_NODE);
  
  //cout << "core_1["<<node1<<"] = " <<core_1[node1] << endl;
  //cout << "core_2["<<node2<<"] = " <<core_2[node2] << endl;

  Node1 attr1 = g1->GetNodeAttr(node1);
  Node2 attr2 = g2->GetNodeAttr(node2);

  if(!nf(attr1[1], attr2[1])){
	  //cout << "here1" << endl;
	  return false;
  }
  
  if(g1->InEdgeCount(node1) > g2->InEdgeCount(node2)
    || g1->OutEdgeCount(node1) > g2->OutEdgeCount(node2)){
	  //cout << "here2" << endl;
    return false;
  }
  
  int i, other1, other2, c_other;
  Edge1 eattr1;
  Edge2 eattr2;
  int termout2=0, termin2=0, new2=0;
  memset(termin2_c,0,classes_count*sizeof(int));
  memset(termout2_c,0,classes_count*sizeof(int));
  memset(new2_c,0,classes_count*sizeof(int));
  

  // Check the 'out' edges of node1
  for(i=0; i<g1->OutEdgeCount(node1); i++)
    { other1=g1->GetOutEdge(node1, i, eattr1);
      c_other = class_1[other1];
      //cout << "other1: " << other1 << endl;
      if (core_1[other1] != NULL_NODE)
        { other2=core_1[other1];
        //cout << "other2: " << other2 << endl;
          if (!g2->HasEdge(node2, other2, eattr2) ||
              !ef(eattr1, eattr2)){
        	  //cout << "here3" << endl;
            return false;
          }
        }
    }
  
  // Check the 'in' edges of node1
  for(i=0; i<g1->InEdgeCount(node1); i++)
    { other1=g1->GetInEdge(node1, i, eattr1);
      c_other = class_1[other1];
      if (core_1[other1]!=NULL_NODE)
        { other2=core_1[other1];
          if (!g2->HasEdge(other2, node2, eattr2) ||
              !ef(eattr1, eattr2)){
        	  //cout << "here4" << endl;
        	  return false;
          }
        }
    }
  
  
  // Check the 'out' edges of node2
  //cout << "checking out edges"  << endl;
  for(i=0; i<g2->OutEdgeCount(node2); i++)
    { other2=g2->GetOutEdge(node2, i);
    //cout << "other2: " << other2 << endl;
      c_other = class_2[other2];
      if (core_2[other2]!=NULL_NODE)
        { 
          other1=core_2[other2];
        	
          //cout << "other1: " << other1 << endl;
          

          /*
           * comment the following lines for checking monomorphism
           *
           */

          if (!g1->HasEdge(node1, other1)){
        	  //cout << "here5" << endl;
        	  return false;
          }

        }
      else
        {
          if (in_2[other2]){
        	  //cout << "in_2["<<other2<<"] is true" <<endl;
          termin2++;
          termin2_c[c_other]++;
          }
          if (out_2[other2]){
        	  //cout << "out_2["<<other2<<"] is true" <<endl;
            termout2++;
            termout2_c[c_other]++;
          }
          if (!in_2[other2] && !out_2[other2]){
            new2++;
            new2_c[c_other]++;
            //cout << "new2_c[" << c_other << "]: " << new2_c[c_other] << endl;
          }
        }
    }
  
  // Check the 'in' edges of node2
  //cout << "checking in edges"  << endl;
  for(i=0; i<g2->InEdgeCount(node2); i++)
    { other2=g2->GetInEdge(node2, i);
    //cout << "other2: " << other2 << endl;
      c_other = class_2[other2];
      if (core_2[other2] != NULL_NODE)
        { 
          other1=core_2[other2];
          //cout << "other1: " << other1 << endl;


          /*
           * comment the following lines for checking monomorphism
           *
           */
          
          if (!g1->HasEdge(other1, node1)){
        	  //cout << "here6" << endl;
        	  return false;
          } 

        }
      else
        { if (in_2[other2]){
        	//cout << "in_2["<<other2<<"] is true" <<endl;
          termin2++;
          termin2_c[c_other]++;
        }
          if (out_2[other2]){
        	  //cout << "out_2["<<other2<<"] is true" <<endl;
            termout2++;
            termout2_c[c_other]++;
          }
          if (!in_2[other2] && !out_2[other2]){
            new2++;
            new2_c[c_other]++;
            //cout << "new2_c[" << c_other << "]: " << new2_c[c_other] << endl;
          }
        }
    }
  
  //Look-ahead check
  if(termin1[core_len] <= termin2 && termout1[core_len] <= termout2){
    for(i = 0; i < classes_count; i++){
      if(termin1_c[core_len][i] > termin2_c[i] ||
         termout1_c[core_len][i] > termout2_c[i]){
    	  //cout << "here7" << endl;
    	  return false;
      }
    }
  }else{
	  //cout << "here8" << endl;
	  return false;
  }
  
  if(new1[core_len] <= new2)
  { //cout << "core_len: " << core_len << endl;
	  for(i = 0; i < classes_count; i++){
		  //cout << "new1_c["<<core_len<<"]["<<i<<"]: " << new1_c[core_len][i] << endl;
		  //cout << "new2_c["<<i<<"]: " << new2_c[i] << endl;
        if(new1_c[core_len][i] > new2_c[i]){
        	//cout << "here9" << endl;
          return false;
        }
      }
  }
  else {
	  //cout << "here10" << endl;
	  return false;
  }
  
  //std::cout << "\nIs Feasible: " << node1 << " " << node2;
  return true;
  
}



/*--------------------------------------------------------------
 * void VF3SubState::AddPair(node1, node2)
 * Adds a pair to the Core set of the state.
 * Precondition: the pair must be feasible
 -------------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
void VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::AddPair(node_id node1, node_id node2)
{
  
  /*std::cout<<"\nAP:";
  print_core(core_1,core_2,n1);
  std::cout<<" <- "<< node1 <<":"<< node2;*/
  
  assert(node1<n1);
  assert(node2<n2);
  assert(core_len<n1);
  assert(core_len<n2);
  assert(class_1[node1] == class_2[node2]);
  
  //Updating the core length
  core_len++;
  added_node1=node1;
  int node_c = class_1[node1];
  core_len_c[node_c]++;
  
  //Checking if node2 is not in T2_in
  if (!in_2[node2])
    { in_2[node2]=core_len;                    // why??
      t2in_len++;
      t2in_len_c[node_c]++;                    // weird: node_c is class of node1
      if (out_2[node2]){
        t2both_len++;
        t2both_len_c[node_c]++;
      }
    }
  
  //Checking if node2 is not in T2_out
  if (!out_2[node2])
    { out_2[node2]=core_len;
      t2out_len++;
      t2out_len_c[node_c]++;
      if (in_2[node2]){
        t2both_len++;
        t2both_len_c[node_c]++;
      }
    }
  
  //Inserting nodes into the core set
  core_1[node1]=node2;
  core_2[node2]=node1;
  
  //Evaluation of the neighborhood
  int i, other, other_c;
  other_c = -1;
  
  for(i=0; i<g2->InEdgeCount(node2); i++)
    { other=g2->GetInEdge(node2, i);
      if (!in_2[other])
        {
        other_c = class_2[other];
        in_2[other]=core_len;
        //in2_set[other_c].push_back(other);
        t2in_len++;
        t2in_len_c[other_c]++;
        if (out_2[other]){
          t2both_len++;
          t2both_len_c[other_c]++;
        }
        }
    }
  
  for(i=0; i<g2->OutEdgeCount(node2); i++)
    { other=g2->GetOutEdge(node2, i);
      if (!out_2[other])
        {
        other_c = class_2[other];
        out_2[other]=core_len;
        //out2_set[other_c].push_back(other);
        t2out_len++;
        t2out_len_c[other_c]++;
        if (in_2[other]){
          t2both_len++;
          t2both_len_c[other_c]++;
        }
        }
    }
  
}



/*--------------------------------------------------------------
 * void VF3SubState::GetCoreSet(c1, c2)
 * Reads the core set of the state into the arrays c1 and c2.
 * The i-th pair of the mapping is (c1[i], c2[i])
 --------------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
void VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::GetCoreSet(node_id c1[], node_id c2[])
{
  int i,j;
  for (i=0,j=0; i<n1; i++)
    if (core_1[i] != NULL_NODE)
      { c1[j]=i;
        c2[j]=core_1[i];
        j++;
      }
}

/*----------------------------------------------------------------
 * Undoes the changes to the shared vectors made by the
 * current state. Assumes that at most one AddPair has been
 * performed.
 ----------------------------------------------------------------*/
template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
void VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::BackTrack()
{
  
  /*std::cout<<"\nBT:";
   print_core(core_1,core_2,n1);
   std::cout<<" -> "<< added_node1 <<":"<< core_1[added_node1];*/
  
  assert(core_len - orig_core_len <= 1);
  assert(added_node1 != NULL_NODE);
  
  int other_c = 0;
  int node_c = class_1[added_node1];
  
  if (orig_core_len < core_len)
    { int i, node2;
      node2 = core_1[added_node1];
      
      if (in_2[node2] == core_len){
        in_2[node2] = 0;
        //in2_set[node_c].erase(node2);
        t2in_len_c[node_c]--;
        if(out_2[node2])
          t2both_len_c[node_c]--;
      }
      
      if (out_2[node2] == core_len){
        out_2[node2] = 0;
        //out2_set[node_c].erase(node2);
        t2out_len_c[node_c]--;
        if(in_2[node2])
          t2both_len_c[node_c]--;
        
      }
      
      //Backtraking neightborhood
      for(i=0; i<g2->InEdgeCount(node2); i++)
        { int other=g2->GetInEdge(node2, i);
          other_c = class_2[other];
          if (in_2[other]==core_len){
            in_2[other]=0;
            //in2_set[other_c].erase(other);
            t2in_len_c[other_c] --;
            if(out_2[other])
              t2both_len_c[other_c]--;
          }
        }
      
      for(i=0; i<g2->OutEdgeCount(node2); i++)
        { int other=g2->GetOutEdge(node2, i);
          other_c = class_2[other];
          if (out_2[other]==core_len){
            out_2[other]=0;
            //out2_set[other_c].erase(other);
            t2out_len_c[other_c] --;
            if(in_2[other])
              t2both_len_c[other_c]--;
          }
        }
      
      core_1[added_node1] = NULL_NODE;
      core_2[node2] = NULL_NODE;
      
      core_len=orig_core_len;
      core_len_c[node_c]--;
      added_node1 = NULL_NODE;
    }
}

template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
bool VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::IsDead() {

  if (n1 > n2){
	  //cout << "first" << endl;
    return true;
  }
  
  if(t1both_len[core_len] > t2both_len ||
     t1out_len[core_len] > t2out_len ||
     t1in_len[core_len] > t2in_len){
	  /*cout << "second" << endl;
	  cout << "core_len: " << core_len << endl;
	  cout <<"t1both_len["<<core_len<<"]: " <<t1both_len[core_len] << ", t2both_len: " << t2both_len << endl;
	  cout <<"t1out_len["<<core_len<<"]: " <<t1out_len[core_len] << ", t2out_len: " << t2out_len << endl;
	  cout <<"t1in_len["<<core_len<<"]: " <<t1in_len[core_len] << ", t2in_len: " << t2in_len << endl;*/
    return true;
  }
  
  for(int c = 0; c < classes_count; c++){
    if(t1both_len_c[core_len][c] > t2both_len_c[c] ||
       t1out_len_c[core_len][c] > t2out_len_c[c] ||
       t1in_len_c[core_len][c] > t2in_len_c[c]){
    /*cout << "third" << endl;
  	  cout << "core_len: " << core_len << endl;
  	  cout << "t1both_len_c[" << core_len << "]["<< c << "]: " << t1both_len_c[core_len][c] << ", t2both_len_c: " << t2both_len_c[c] << endl;
  	  cout << "t1out_len_c[" << core_len << "]["<< c << "]: " << t1out_len_c[core_len][c] << ", t2out_len_c: " << t2out_len_c[c] << endl;
  	  cout << "t1in_len_c[" << core_len << "]["<< c << "]: " << t1in_len_c[core_len][c] << ", t2in_len_c: " << t2in_len_c[c] << endl;*/
      return true;
    }
  }
  
  return false;
}


template <typename Node1, typename Node2,
typename Edge1, typename Edge2,
typename NodeComparisonFunctor, typename EdgeComparisonFunctor>
VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>* VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>::GetParent()
{
  return (VF3SubState<Node1,Node2,Edge1,Edge2,NodeComparisonFunctor,EdgeComparisonFunctor>*)parent;
}


#endif
