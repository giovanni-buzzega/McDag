#include <iostream>
#include <sys/time.h> // measure time
#include <sstream>	// for easy string construction
#include <vector>
#include <cassert>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <math.h>
#include <random>
#include <map>
typedef  int node_id_t ; // using -1 as a value, do not unsign 
using namespace std;


const bool DEBUG = false;
bool d_flag = false;
timeval t_begin, t_end;

double time_elapsed(const timeval begin, const timeval end) {
	  return end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1e-6;
}

class Match : public vector<int> {
public:
	node_id_t pos;
	Match(const std::vector<int>& v,int assigned_pos) : vector<int>(v) {
		pos = assigned_pos;
	}
	Match(const std::vector<int>& v) : vector<int>(v) {
		pos = -1;
	}

	const std::string to_string() const {
		ostringstream os;
		os <<"(";
		for (size_t i = 0; i < this->size()-1; i++)
			os << this->at(i) << ", ";
		os << this->back() << ")";
		return os.str();
	}
	friend ostream &operator<<(ostream &os, const Match &m) {
		os << m.to_string();
		return os;
	}
	bool operator<(const Match& other) const {
		assert(this->size() == other.size());
		for (size_t i = 0; i < this->size(); i++)
			if (!(this->at(i) < other.at(i)))
				return false;
		return true;
	}
	// domination O(K)
	bool operator>(const Match& other) const {
		assert(this->size() == other.size());
		for (size_t i = 0; i < this->size(); i++)
			if (!(this->at(i) > other.at(i)))
				return false;
		return true;
	}
	bool operator<=(const Match& other) const {
		assert(this->size() == other.size());
		for (size_t i = 0; i < this->size(); i++)
			if (this->at(i) > other.at(i))
				return false;
		return true;
	}
	bool operator>=(const Match& other) const {
		assert(this->size() == other.size());
		for (size_t i = 0; i < this->size(); i++)
			if (this->at(i) < other.at(i))
				return false;
		return true;
	}
	bool operator==(const Match& other) const {
		assert(this->size() == other.size());
		for (size_t i = 0; i < this->size(); i++)
			if (this->at(i) != other.at(i))
				return false;
		return true;
	}
};
// each match is wide K; matches is 
vector<int> fishlet(const vector<Match*> & matches) {
	assert(matches.size() > 0);
	vector<int> final_result(*(matches[0]));
	for (size_t i = 1; i < matches.size(); i++) 
		for (size_t j = 0; j < matches[i]->size(); j++) 
			if (final_result[j] > matches[i]->at(j))
				final_result[j] = matches[i]->at(j);
	return final_result; 
}

class Final_match : public Match{
public:
	Final_match(const vector<Match*>& m) : Match(fishlet(m)) {
		modeled.reserve(m.size()); 
		for (size_t i = 0; i < m.size(); i++)
			modeled.push_back(m[i]->pos);
		pos=-1;
	}
	// This should always be sorted 
	vector<node_id_t> modeled;
	bool operator==(const Final_match& other) const {
		if (this->modeled.size() != other.modeled.size()) return false;
		for (size_t i = 0; i < this->modeled.size(); i++)
			if (this->modeled[i] != other.modeled[i])
				return false;
		return true;
	}
};
// (todo -> change to universal hashing)
namespace std {
// Hash for Match; sax_hash (shift add xor)
template <> struct hash<Match> {
	size_t operator()(const Match & m) const {
		size_t h = 0;

		for (int i : m) {
			h ^= (h << 5) + (h >> 2) + i;
		}
		return h;
	}
};


// Hash for minimization
class Hash_node {
private:
	constexpr uint64_t universal_hash_64(const uint64_t x) {    // hashes x universally into 32 + 32 bits using random ODD seed a1, a2 > 0
		//constexpr uint64_t a1 = 0x4f2162926e40c299UL;
		//constexpr uint64_t a2 = 15429903764174996207UL;
		//return (static_cast<uint64_t>((x ^ (x >> 33U)) * a1) & 0xffffffff00000000UL) | (static_cast<uint64_t>((x ^ (x >> 33U)) * a2) >> 32U);    // see https://en.wikipedia.org/wiki/Universal_hashing
		return (static_cast<uint64_t>((x ^ (x >> 33U)) * 0x4f2162926e40c299UL) & 0xffffffff00000000UL) | (static_cast<uint64_t>((x ^ (x >> 33U)) * 15429903764174996207UL) >> 32U);    // see https://en.wikipedia.org/wiki/Universal_hashing
	}

public:
	uint64_t operator()(const vector<uint64_t> & x) {
		uint64_t res = 0;
		for (size_t i = 0; i < x.size(); i++){
			res = res ^ universal_hash_64(x[i]);
		}
		return res;
	}
};

}
namespace std {
// permutation invariant hash, taken from python frozen set https://stackoverflow.com/questions/20832279/python-frozenset-hashing-algorithm-implementation
template <> struct hash<Final_match> {
	size_t operator()(const Final_match & m) const {
		const vector<node_id_t>& vec = m.modeled;
		size_t hash = 1927868237UL;

		hash *= vec.size() + 1;
		for (const node_id_t pos : vec) {
			/* Work to increase the bit dispersion for closely spaced hash
			   values.  The is important because some use cases have many
			   combinations of a small number of elements with nearby
			   hashes so that many distinct combinations collapse to only
			   a handful of distinct hash values. */
			hash ^= (pos ^ (pos << 16) ^ 89869747UL)  * 3644798167UL;
		}
		hash = hash * 69069U + 907133923UL;
		// if (hash == -1)
		//	 hash = 590923713UL;
		return hash;
	}
};
}

template <typename Match> class ST_graph {
private:
	//in d1 I need to index nodes by their coordinates
	//in d2 I need to index nodes by their modeled nodes
	long double num_paths; 
	vector < long double> sum_out_deg; // this is for recursively counting paths
	vector<long double> hist_lengths;
	vector<vector<node_id_t>> adj_list;	
	vector<vector<node_id_t>> reverse_adj_list;
	bool deterministic;
	bool minimal;
	bool use_reverse;
	bool use_adj_list;
	size_t num_edges;
	vector<Match> node_list;
	size_t num_nodes;
	unordered_map<Match , node_id_t> node_lookup; 
	// Building: take a node: look at the modeled (vector of ints or of match*?),  lookup children of modeled (better int) , group them (better match*, for fishlet) (n*sigma for domination filter), build new bitmap, add edge/node
public:
	ST_graph( bool adj, bool reverse, bool det) { 
		minimal = false;
		num_edges = 0;
		num_nodes = 0;
		num_paths = 0;
		use_adj_list = adj; 
		use_reverse = reverse; 
		deterministic = det;
		assert ((deterministic && use_reverse) || (!deterministic && use_adj_list));
		}; // TODO reserve UB spots
	node_id_t source;
	node_id_t sink;

	void dfs(const node_id_t node, const string & prefix,const char final_letter, const vector<string> & strings){
		char mychar = strings[0][get_node(node)[0]];
		string newprefix = prefix+mychar;
		if (mychar == final_letter)
			cout << newprefix <<endl;
		for (node_id_t child: get_children_of(node))
			dfs(child, newprefix, final_letter, strings);
	}

	// TODO: this can be computed during construction without using reverse_adj_list (with priority queue)
	void ratio_path_lengths(const vector<string> & strings, long double & card_lcs, long double & card_lcs_1, int & lcs_length) {
		size_t n = strings[0].size(); // n is the minimum length of the string (aka maximum length of any cs)
		for (const string & s : strings)
			if (s.size() < n)
				n = s.size();
		if (!hist_lengths.empty())  {
			lcs_length = hist_lengths.size();
			card_lcs = hist_lengths[lcs_length-1];
			card_lcs_1 = hist_lengths[lcs_length-2];
			return;
		}
		vector<long double> num_lcs(adj_list.size(),0);
		vector<long double> num_lcs_1(adj_list.size(),0);
		vector<int> len_lcs(adj_list.size(),0);
		// Add a node when all children have been visited
		vector<int> children_left; 
		children_left.resize(adj_list.size());
		for (size_t i = 0; i < adj_list.size(); i++) 
			children_left[i] = adj_list[i].size();
		// Queue of nodes that have 0 children_left
		queue<node_id_t> q; 

		q.push(sink);
		len_lcs[sink] = 1;
		num_lcs[sink] = 1;

		for (; !q.empty(); q.pop()) {
			node_id_t curr_node = q.front();
			// Update based on the value of the children
			for (node_id_t child: adj_list.at(curr_node)) {
				if ( len_lcs[curr_node] == len_lcs[child]) {
					num_lcs_1[curr_node] = num_lcs_1[child]+num_lcs[curr_node];
					num_lcs[curr_node] = num_lcs[child];
					len_lcs[curr_node] = len_lcs[child]+1;
				} else if ( len_lcs[curr_node] < len_lcs[child]) {
					num_lcs_1[curr_node] = num_lcs_1[child];
					num_lcs[curr_node] = num_lcs[child];
					len_lcs[curr_node] = len_lcs[child]+1;
				} else if ( len_lcs[curr_node] == len_lcs[child]+1) {
					num_lcs_1[curr_node] += num_lcs_1[child];
					num_lcs[curr_node] += num_lcs[child];
				} else if ( len_lcs[curr_node] == len_lcs[child]+2) {
					num_lcs_1[curr_node] += num_lcs[child];
				}
			}
			// Add parent if visiting its last children
			for (node_id_t parent: reverse_adj_list.at(curr_node)) {
				children_left[parent]--;
				if (children_left[parent] == 0) {
					q.push(parent);
				}
			}
		}

		card_lcs = num_lcs[source];
		card_lcs_1 = num_lcs_1[source];
		lcs_length = len_lcs[source]-2;
		return;
	}
	// O(largest_antichain * n^2)
	void histogram_path_lengths() {
		if (!hist_lengths.empty()) return ;
		
		vector<vector<long double>> hist_of_node;
		// when all parents have been visited, clear memory of hist_of_nodes
		vector<int> parents_left;
		// Add a node when all children have been visited
		vector<int> children_left; 
		parents_left.resize(adj_list.size());
		children_left.resize(adj_list.size());
		hist_of_node.resize(adj_list.size());
		for (size_t i = 0; i < adj_list.size(); i++) 
			children_left[i] = adj_list[i].size();
		// Queue of nodes that have 0 children_left
		queue<node_id_t> q; 

		hist_of_node.insert(hist_of_node.cbegin() + sink,vector<long double>(1,(long double) 1));
		//TODO: reserve for each hist of node;
		q.push(sink);

		for (; !q.empty(); q.pop()) {
			node_id_t curr_node = q.front();
			parents_left[curr_node] = reverse_adj_list.at(curr_node).size();
			// Update based on the value of the children
			for (node_id_t child: adj_list.at(curr_node)) {
				for (size_t i = 0; i < hist_of_node[child].size(); i++) 
					hist_of_node[curr_node][i+1] += hist_of_node[child][i];
				parents_left[child]--; 
				// clear memory
				if (parents_left[child]==0) {
					hist_of_node[child].clear();
				}
			}
			// Add parent if visiting its last children
			for (node_id_t parent: reverse_adj_list.at(curr_node)) {
				children_left[parent]--;
				if (children_left[parent] == 0) {
					size_t max_child_length = hist_of_node[curr_node].size();
					for (const node_id_t child: adj_list.at(parent)) 
						if (hist_of_node[child].size() > max_child_length)
							max_child_length = hist_of_node[child].size();
					hist_of_node[parent] = vector<long double>(max_child_length+1,0);
					q.push(parent);
				}
			}
		}

		hist_of_node[source].erase(hist_of_node[source].begin());
		hist_lengths = hist_of_node[source];
	}
	// O(largest_antichain * n)
	long double iterative_count_paths() {
		if (num_paths > 0) return num_paths;
		if (!hist_lengths.empty()) { 
			long double c = 0;
			for (size_t i = 0; i < hist_lengths.size(); i++)
				c+= hist_lengths[i];
			return c;
		}

		vector< long double> sum_out_deg; // this is for counting paths
		// when all parents have been visited, clear memory of sum_out_deg
		// Add a node when all children have been visited
		vector<int> children_left; 
		children_left.resize(adj_list.size());
		sum_out_deg.resize(adj_list.size());
		for (size_t i = 0; i < adj_list.size(); i++) 
			children_left[i] = adj_list[i].size();
		// Queue of nodes that have 0 children_left
		queue<node_id_t> q; 

		sum_out_deg[sink] = 1;
		q.push(sink);

		for (; !q.empty(); q.pop()) {
			node_id_t curr_node = q.front();
			// Update based on the value of the children
			for (node_id_t child: adj_list.at(curr_node)) {
				sum_out_deg[curr_node] += sum_out_deg[child];
			}
			// Add parent if visiting its last children
			for (const node_id_t parent: reverse_adj_list.at(curr_node)) {
				children_left[parent]--;
				if (children_left[parent] == 0) {
					sum_out_deg[parent] = 0;
					q.push(parent);
				}
			}
		}
		num_paths = sum_out_deg[source];
		return num_paths;
	}
	size_t get_num_nodes() { 
		if (num_nodes == 0)
			num_nodes = node_list.size();
		return num_nodes;
	}
	size_t get_num_edges() { 
		if (num_edges == 0) {
			if (use_adj_list) 
				for (size_t i = 0; i < adj_list.size(); i++)
					num_edges+= adj_list[i].size();
			else for (size_t i = 0; i < reverse_adj_list.size(); i++)
					num_edges+= reverse_adj_list[i].size();
		}
		return num_edges; }
	Match & get_source() { return node_list[source]; }
	Match & get_sink() { return node_list[sink]; }
	Match & get_node(const node_id_t p) { return node_list[p]; }
	vector<node_id_t> & get_children_of(const Match * match) { assert(use_adj_list); return adj_list[match->id]; }
	vector<node_id_t> & get_children_of(const int id) { assert(use_adj_list); return adj_list[id]; }
	vector<node_id_t> & get_parents_of(const int id) { assert(use_reverse); return reverse_adj_list[id]; }
	bool contains_match(const Match & m) { 
		if (m.pos != -1) return true; // if position is populated it should be in the graph
		else return !(node_lookup.find(m) == node_lookup.end());}
	node_id_t add_node(const Match & n) {
		assert (!contains_match(n));
		int new_pos = node_list.size();
		node_list.push_back(n); 
		num_nodes++;
		node_list[new_pos].pos = new_pos;
		if (use_adj_list)
			adj_list.push_back(vector<node_id_t>());
		if (use_reverse)
			reverse_adj_list.push_back(vector<node_id_t>());
		node_lookup.insert({n,new_pos});
		return new_pos;
	} 
	node_id_t get_pos(const Match & n) {
		assert (contains_match(n));
		return node_lookup.find(n)->second;
	} 
	void add_edge(const node_id_t match, const node_id_t child) {
		if (use_adj_list)
			adj_list[match].push_back(child);
		if (use_reverse)
			reverse_adj_list[child].push_back(match);
	}
	void add_edge(const Match & match, const Match & child) {
		if (use_adj_list)
			adj_list[match.pos].push_back(child.pos);
		if (use_reverse)
			reverse_adj_list[child.pos].push_back(match.pos);
	}
int collision_detector(vector<uint64_t> & hashes , vector<vector<node_id_t>> & prev_adj, node_id_t curr_node, node_id_t prev_node, char c_curr, char c_prev, bool virt) {
	// TODO: remove virtual minimization as it is not a bottleneck
	//if virt we don't update all edges, so the graph is really consistent
		// Collision here must be verified via hashes
	if (c_curr != c_prev) {
		cout << "COLLISION DETECTED (different letter) " << prev_node << get_node(prev_node)
			<< " <- " << curr_node << get_node(curr_node) << " # " << hashes[prev_node] << " # " << hashes[curr_node] << endl;
		return -2;
	}

	if (prev_adj[curr_node].size() != prev_adj[prev_node].size()) {
		cout << "DIFFERENT SIZES "<<prev_adj[curr_node].size() << " " <<prev_adj[prev_node].size() <<endl;
		return -3;
	}
	for(node_id_t prev_neigh_c : prev_adj[curr_node]) {
		bool found = false;
		for(node_id_t pnn : prev_adj[prev_node]) {
			if ((virt && hashes[pnn] == hashes[prev_neigh_c]) || 
					pnn == prev_neigh_c) {
				found = true;
				break;
			}
		}
		if (!found) {
			cout << "COLLISION DETECTED (child not found in the right node) " << prev_node 
				<< get_node(prev_node) << " <- " << curr_node << get_node(curr_node) 
				<< " nf: "<< prev_neigh_c << get_node(prev_neigh_c) << " #: "<< hashes[prev_neigh_c] << endl;
			cout << "hashes right node children:";
			for (node_id_t pnc: prev_adj[curr_node])
				cout << " # " << hashes[pnc] << " nf " << pnc << get_node(pnc);
			cout <<endl;
			return -1;
		}
	}
	return 0;
}
	bool minimalize(const vector<string>& strings, const vector<char> & sigma, bool virt) {
		assert (use_adj_list && use_reverse);
		if (virt && minimal)
			return true;
		get_num_nodes();
		get_num_edges();
		// Revuz's algorithm for dag minimization
		// queue of non-processed nodes whose out_neighbors have been processed (in neighbors for the non-deterministic case)
		vector<node_id_t> frontier;
		vector<node_id_t> next_frontier;
		vector<uint64_t> hashes;

		vector<vector<node_id_t>> & next_adj = deterministic? reverse_adj_list: adj_list;
		vector<vector<node_id_t>> & prev_adj = !deterministic? reverse_adj_list: adj_list;

		vector<int> neighbors_left; 
		neighbors_left.resize(node_list.size());
		hashes.resize(node_list.size());
		for (size_t i = 0; i < node_list.size(); i++)  
			neighbors_left[i] = prev_adj[i].size();

		Hash_node hasher;
		if (deterministic) {// start from sink 
			next_frontier.push_back(sink);
			hashes[sink] =	hasher(vector<uint64_t>(sigma.size()+1,(uint64_t)'$'));
		} else { // start from source
			next_frontier.push_back(source);
			hashes[source] = hasher(vector<uint64_t>(sigma.size()+1,(uint64_t)'#'));
		}

		int h = 0; //height defined as the max distance from the sink
		while(!next_frontier.empty()) {
			// cout << "height "<<h<<endl;
			h++;
			frontier = std::move(next_frontier);
			for (size_t i = 0; i < frontier.size(); i++) {
				node_id_t curr_node = frontier[i];
				// compute hash
				vector<uint64_t> neigh_hashes(sigma.size()+1,0);
				for(node_id_t prev_neigh : prev_adj[curr_node]) 
					neigh_hashes.push_back(hashes[prev_neigh]);
				// also hash my character
				neigh_hashes.push_back((uint64_t)strings[0][get_node(curr_node)[0]]);
				hashes[curr_node] = hasher(neigh_hashes);
				
				// Add neighbor if you are its last visited prev_neighbor
				for (node_id_t neigh: next_adj.at(curr_node)) {
					neighbors_left[neigh]--;
					if (neighbors_left[neigh] == 0) 
						next_frontier.push_back(neigh);
				}
			}
			// unordered_map<uint64_t, vector<node_id_t>> hash_map;
			// // map<uint64_t, vector<node_id_t>> hash_map;
			// for (size_t i = 0; i < frontier.size(); i++) {
			// 	node_id_t curr_node = frontier[i];
			// 	hash_map[hashes[curr_node]].push_back(curr_node);
			// 	if (hash_map[hashes[curr_node]].size() > 1) {
			// 		node_id_t prev_node = hash_map[hashes[curr_node]][hash_map[hashes[curr_node]].size()-2]; 
			// 		char c_curr = strings[0][get_node(curr_node)[0]];
			// 		char c_prev = strings[0][get_node(prev_node)[0]];
			// 		if (collision_detector(hashes, prev_adj, curr_node, prev_node,c_curr, c_prev,virt )!= 0)
			// 			return false;
			// 	}
			// }
			// merge(hash_map,virt);
			
			// sort signature
			auto my_comp = [&hashes](const node_id_t i1, const node_id_t i2)
			     { return hashes[i1] < hashes[i2];};
			sort(frontier.begin(),frontier.end(),my_comp); // n log n
			// merge equal nodes (they must have the same letter and point to the same nodes)
			for (size_t i = frontier.size()-1; i > 0; i--) {
				const node_id_t curr_node = frontier[i];
				const node_id_t prev_node = frontier[i-1];
				if (hashes[curr_node] == hashes[prev_node]) {
					// i-1 <- i
					char c_curr = strings[0][get_node(curr_node)[0]];
					char c_prev = strings[0][get_node(prev_node)[0]];
					if (collision_detector(hashes, prev_adj, curr_node, prev_node,c_curr, c_prev,virt )!= 0)
						return false;
					merge(prev_node,curr_node,virt);
				}
			}
		}
		minimal = true;
		return true;
	}

	void merge( const unordered_map<uint64_t, vector<node_id_t>> & hash_map, const bool virt) {
		vector<vector<node_id_t>> & det_adj = deterministic?  adj_list : reverse_adj_list;
		vector<vector<node_id_t>> & non_det_adj = !deterministic? adj_list : reverse_adj_list;
		for (auto const& pair : hash_map) {
			node_id_t m1 = pair.second[0];

			// parents of m2 just change child, so no edge is lost
			// children of m2 lose an edge if they are children of m1 as well
			// Since we are merging this is guaranteed
			num_edges -= (det_adj[m1].size() * (pair.second.size()-1));
			num_nodes -= (pair.second.size()-1);
			if (!virt) {
				for (const node_id_t ch_pos : det_adj[m1]) {
					// update children to remove all nodes to be merged
					for (int i = non_det_adj[ch_pos].size()-1; i >= 0; i--) {
						if (non_det_adj[ch_pos][i] != m1 && 
								find(pair.second.begin(), pair.second.end(),non_det_adj[ch_pos][i]) != pair.second.end())
							non_det_adj[ch_pos].erase(non_det_adj[ch_pos].begin()+i);
					}
				}
			}
			for (size_t i = 1; i < pair.second.size(); i++) {
				node_id_t m2 = pair.second[i];
				if (!virt) {
					for (const node_id_t par_pos : non_det_adj[m2]) {
						// update parents to take new child
						*find(det_adj[par_pos].begin(), det_adj[par_pos].end(),m2) = m1;
						// update new child to take parents
						non_det_adj[m1].push_back(par_pos);
					}
				}
				// delete link to parents
				non_det_adj[m2].clear();
				// delete link to children
				det_adj[m2].clear();
			}
		}
	}
	void merge(const node_id_t m1, const node_id_t m2, bool virt) {
		// m1 <- m2
		vector<vector<node_id_t>> & det_adj = deterministic?  adj_list : reverse_adj_list;
		vector<vector<node_id_t>> & non_det_adj = !deterministic? adj_list : reverse_adj_list;
		
		// parents of m2 just change child, so no edge is lost
		// children of m2 lose an edge if they are children of m1 as well
		// Since we are merging this is guaranteed for all children
		num_edges -= det_adj[m2].size();
		num_nodes--;
		if (!virt) {
			for (const node_id_t par_pos : non_det_adj[m2]) {
				*find(det_adj[par_pos].begin(), det_adj[par_pos].end(),m2) = m1;
				non_det_adj[m1].push_back(par_pos);
			}
			for (const node_id_t ch_pos : det_adj[m2]) {
				non_det_adj[ch_pos].erase(find(non_det_adj[ch_pos].begin(), non_det_adj[ch_pos].end(),m2));
			}
		}
		// delete link to parents
		non_det_adj[m2].clear();
		// delete link to children
		det_adj[m2].clear();
	}
	int print_stats(bool z_flag, bool m_flag, bool l_flag, bool r_flag, bool virt, const string & name, const vector<string>& strings, const vector<char> & sigma, bool do_intensive_computations) {
		cout << "number of "<<name<<" paths: "<< iterative_count_paths() <<endl;
		if (z_flag && do_intensive_computations) {
			if (!minimalize(strings,sigma,virt)) {
				cout << "minimization of "<<name<<" failed. exit\n";
				return -1;
			}
			gettimeofday(&t_end, 0);
			cout << "time for d1r: " << time_elapsed(t_begin, t_end) << " seconds."<< endl;
			cout << "number of "<<name<<"_minimal nodes: "<< get_num_nodes() << ", number of "<<name<<"_minimal edges: "<< get_num_edges()<<endl;
		}
		if (m_flag && do_intensive_computations) {
			cout<< "language of "<<name<<":\n";
			dfs(source, "", '$', strings);
		}
		if (l_flag && do_intensive_computations) {
			cout<< "Distribution of length of "<<name<<" paths:\n";
			histogram_path_lengths();
			for (size_t i = 0; i < hist_lengths.size(); i++)
				if (hist_lengths[i] > 0)
					cout<<i<<": "<<hist_lengths[i]<<endl;
		}
		if (r_flag) {
			long double card_lcs = 0;
			long double card_lcs_1 = 0;
			int lcs_length = 0;
			ratio_path_lengths( strings, card_lcs, card_lcs_1, lcs_length) ;
			cout<< "Ratio according to "<<name<<" #(lcs-1): " << card_lcs_1 << " #lcs: "<<card_lcs << " lcs_length: "<<lcs_length <<endl;
		}
		return 0;
	}
};


// O(|sigma|^2) in general (for 2 strings it can be |sigma|log|sigma|) 
vector<Match> next_matches(const Match & match, const vector<string>& strings, unordered_map<char,vector<vector<int>>>& occurrence, const unordered_set<char>& allowed_chars, bool rightmost) {
	const int K = strings.size();
	// Create one match per letter
	vector<Match> candidate_matches;
	candidate_matches.reserve(allowed_chars.size());
	for (char c: allowed_chars) {
		vector<int> next_match(K);
		int i;
		for (i = 0; i < K; i++) {
			next_match[i] = occurrence[c][i][match[i]];
			if (next_match[i] == -1) i = K+1;	// do not create any prev match
		}
		if (i <= K) candidate_matches.push_back(Match(next_match));
	}
	// Filter non-rightmost (non-leftmost) ones
	if (K != 2 || candidate_matches.size() < 10) {
		for (size_t i = 0; i < candidate_matches.size(); i++) {
			for (size_t j = i+1; j < candidate_matches.size(); j++) {
				if ((rightmost && candidate_matches[i] < candidate_matches[j]) ||
						(!rightmost && candidate_matches[i] > candidate_matches[j])) {
					if (DEBUG) cout << "removing "<< i << " " << candidate_matches[i][0] << " " << candidate_matches[i][1] <<endl;
					if (DEBUG) cout << "because of "<< j << " " << candidate_matches[j][0] << " " << candidate_matches[j][1] <<endl;
					candidate_matches.erase(candidate_matches.begin()+i);
					i--; j=candidate_matches.size();
				}
				else if ((rightmost && candidate_matches[i] > candidate_matches[j] ) ||
						(!rightmost && candidate_matches[i] < candidate_matches[j] )){	 
					if (DEBUG) cout << "removing "<< j << " " << candidate_matches[j][0] << " " << candidate_matches[j][1] <<endl;
					if (DEBUG) cout << "because of "<< i << " " << candidate_matches[i][0] << " " << candidate_matches[i][1] <<endl;
					candidate_matches.erase(candidate_matches.begin()+j);
					j--;
				}
			}
		}
	} else { // more efficient filter for K=2 // more efficient in theory, not in practice as sort has a huge overhead for low sigma values
		auto my_comp = (rightmost) ? 
			[](const Match & l, const Match & r)
			{ return  (l[0] > r[0] || (l[0] == r[0] && l[1] > r[1]));}:
			[](const Match & l, const Match & r)
			{ return  (l[0] < r[0] || (l[0] == r[0] && l[1] < r[1]));};
		sort(candidate_matches.begin(), candidate_matches.end(), my_comp);
		int best_so_far = candidate_matches[0][1];
		for (size_t i = 1; i < candidate_matches.size(); i++) {
			if ((rightmost && candidate_matches[i][1] < best_so_far) || 
					(!rightmost && candidate_matches[i][1] > best_so_far)) {
				if (DEBUG) cout<<"erase "<< i<< " " << candidate_matches[i]<<endl;
				candidate_matches.erase(candidate_matches.begin()+i);
				i--;
			} else
				best_so_far = candidate_matches[i][1];
		}
	}
	return candidate_matches;
}
void build_csa(ST_graph<Match> *csa_ptr,const vector<string> & strings, unordered_map<char,vector<vector<int>>> & occurrence, const vector<char> & sigma, bool codeterministic) {
	const unordered_set<char> sigma_set(sigma.begin(),sigma.end());
	int K = strings.size();
	ST_graph<Match> & csa = *csa_ptr;
	// Create D1 and add its sink
	vector<int> m(K);
	const char c = strings[0].back();
	for (int i = 0; i < K; i++) {
		assert(strings[i].back() == c);
		m[i] = strings[i].size()-1;
	}
	csa.source = csa.add_node(Match(vector<int>(K,0)));
	csa.sink = csa.add_node(Match(m));

	queue<node_id_t> q;
	if (codeterministic) q.push(csa.sink);
	else q.push(csa.source);

	for (; !q.empty(); q.pop()) {
		node_id_t curr_node = q.front();
		if (DEBUG) cout<<curr_node<<csa.get_node(curr_node).to_string()<<endl;

		Match & curr_match = csa.get_node(curr_node);
		// Create one match per letter
		vector<Match> next;
		next.reserve(sigma.size());
		for (char c: sigma) {
			vector<int> next_match(K);
			int i;
			for (i = 0; i < K; i++) {
				next_match[i] = occurrence[c][i][curr_match[i]];
				if (next_match[i] == -1) i = K+1;	// do not create any prev match
			}
			if (i <= K) next.push_back(Match(next_match));
		}


		if (codeterministic) csa.add_edge(csa.source,curr_node);
		else  csa.add_edge(curr_node,csa.sink);
		if (DEBUG) {cout << "adding edge "<<csa.get_source().to_string() << "->" << csa.get_node(curr_node).to_string() <<"\n";}
		for (Match & m: next) {
			if (!(csa.contains_match(m))) {
				q.push(csa.add_node(m));
			}
			if (DEBUG) cout << "adding edge "<<m.to_string() << "->" << csa.get_node(curr_node).to_string() <<"\n";
			if (codeterministic) csa.add_edge(csa.get_pos(m),curr_node);
			else  csa.add_edge(curr_node,csa.get_pos(m));
		}
	}
}
void build_d1(ST_graph<Match> *d1_ptr,const vector<string> & strings, unordered_map<char,vector<vector<int>>> & occurrence, const vector<char> & sigma, bool codeterministic) {
	const unordered_set<char> sigma_set(sigma.begin(),sigma.end());
	int K = strings.size();
	ST_graph<Match> & d1 = *d1_ptr;
	// Create D1 and add its sink
	vector<int> m(K);
	const char c = strings[0].back();
	for (int i = 0; i < K; i++) {
		assert(strings[i].back() == c);
		m[i] = strings[i].size()-1;
	}
	d1.source = d1.add_node(Match(vector<int>(K,0)));
	d1.sink = d1.add_node(Match(m));

	queue<node_id_t> q;
	if (codeterministic) q.push(d1.sink);
	else q.push(d1.source);

	for (; !q.empty(); q.pop()) {
		node_id_t curr_node = q.front();
		if (DEBUG) cout<<curr_node<<d1.get_node(curr_node).to_string()<<endl;
		vector<Match> next = next_matches(d1.get_node(curr_node), strings, occurrence, sigma_set, codeterministic); 
		if (next.size() == 0) {
			if (codeterministic) d1.add_edge(d1.source,curr_node);
			else  d1.add_edge(curr_node,d1.sink);
			if (DEBUG) cout << "adding edge "<<d1.get_source().to_string() << "->" << d1.get_node(curr_node).to_string() <<"\n";}
		for (Match & m: next) {
			if (!(d1.contains_match(m))) {
				q.push(d1.add_node(m));
			}
			if (DEBUG) cout << "adding edge "<<m.to_string() << "->" << d1.get_node(curr_node).to_string() <<"\n";
			if (codeterministic) d1.add_edge(d1.get_pos(m),curr_node);
			else  d1.add_edge(curr_node,d1.get_pos(m));
		}
	}
}

void codeterminize_d1(ST_graph<Match> *d1_codet_ptr, ST_graph<Match> * d1_ptr,const vector<string> & strings,unordered_map<char,vector<vector<int>>> & occurrence, const vector<char> & sigma) {
	const unordered_set<char> sigma_set(sigma.begin(),sigma.end());
	// Create D1_codet and add its sink
	ST_graph<Match> & d1 = *d1_ptr;
	ST_graph<Match> & d1_codet = *d1_codet_ptr;
	d1_codet.sink = d1_codet.add_node(Match(vector<int>(d1.get_sink().cbegin(),d1.get_sink().cend())));
	d1_codet.source = d1_codet.add_node(Match(vector<int>(d1_codet.get_sink().size(),0)));
	// keep track of where we are in d1 (not sure about size d1_codet, so use a map)
	vector< vector<node_id_t>> corresponding_nodes;
	corresponding_nodes.resize(d1.get_num_nodes());

	auto my_comp = [&d1_codet](const node_id_t i1, const node_id_t i2)
	{ Match & m1 = d1_codet.get_node(i1);
	  Match & m2 = d1_codet.get_node(i2);
	  return  (*max_element(m2.cbegin(), m2.cend())) > (*max_element(m1.cbegin(), m1.cend())); };
	priority_queue<node_id_t,
		vector<node_id_t>,
		decltype(my_comp)> q{my_comp};
	q.push(d1_codet.sink);
	corresponding_nodes[d1_codet.sink].push_back(d1.sink);

	for (; !q.empty(); ) {
		node_id_t curr_node = q.top();
		q.pop();
		//find allowed characters based on parents
		unordered_set<char> allowed_chars;
		unordered_map<char, vector<node_id_t>> parents; // by determinism you cannot add the same parent twice
		// cout << "parents for each letter: ";
		for (const node_id_t m: corresponding_nodes[curr_node]) {
			for (const node_id_t par_pos: d1.get_parents_of(m)) {
				Match & parent = d1.get_node(par_pos);
				parents[strings[0][parent[0]]].push_back(par_pos); 
			}
		}
		vector<Match> next = next_matches(d1_codet.get_node(curr_node), strings, occurrence, sigma_set, true);
		// smart filter: filter strictly dominated parents by some rightmost match
		for (auto par_iter = parents.begin(); par_iter != parents.end(); ++par_iter) {
			for (size_t i = 0; i < par_iter->second.size(); i++) {
				Match & parent_1 = d1.get_node(par_iter->second[i]);
				for (size_t j = 0 ; j < next.size(); j++) {
					Match & parent_2 = next[j];
					if (parent_1 < parent_2) {
						par_iter->second.erase(par_iter->second.begin()+i); // TODO: do not use erase
						i--;
						j = next.size();
					}
				}
			}
		}
		// Only parents that survived the filter define allowed characters
		for (auto par_iter = parents.begin(); par_iter != parents.end(); ++par_iter) {
			for (size_t i = 0; i < par_iter->second.size(); i++) {
				Match & parent_1 = d1.get_node(par_iter->second[i]);
				if (strings[0][parent_1[0]] != '#')
					allowed_chars.insert(strings[0][parent_1[0]]);
			}
		}
		next = next_matches(d1_codet.get_node(curr_node), strings, occurrence, allowed_chars, true);

		if (next.size() == 0) {
			d1_codet.add_edge(d1_codet.source,curr_node);
			if (DEBUG) cout << "adding edge "<<d1_codet.get_source().to_string() << "->" << d1_codet.get_node(curr_node).to_string() <<"\n";}
		for (Match & m: next) {
			if (!(d1_codet.contains_match(m))) {
				q.push(d1_codet.add_node(m));
			}
			char curr_char = strings[0][m[0]];
			node_id_t m_pos = d1_codet.get_pos(m);
			if (corresponding_nodes.size() <= (size_t)m_pos)
				corresponding_nodes.resize((size_t)corresponding_nodes.size()*1.2 + 1);
			if (!corresponding_nodes[m_pos].empty() ) {
				int old_size = corresponding_nodes[m_pos].size();
				// do not add parent if it is already present
				for (const node_id_t i : parents[curr_char])
					//if (!std::binary_search(corresponding_nodes[m_pos].begin(), corresponding_nodes[m_pos].begin()+old_size,i))
					if (find(corresponding_nodes[m_pos].begin(), corresponding_nodes[m_pos].begin()+old_size,i)==corresponding_nodes[m_pos].begin()+old_size)
						corresponding_nodes[m_pos].push_back(i);
			} else
				corresponding_nodes[m_pos] = parents[curr_char];
			//std::sort(corresponding_nodes[m_pos].begin(),corresponding_nodes[m_pos].end());
			d1_codet.add_edge(m_pos,curr_node);
			if (DEBUG) cout << "adding edge "<<m.to_string() << "->" << d1_codet.get_node(curr_node).to_string() <<"\n";
		}
		// release memory
		corresponding_nodes[curr_node].clear();
	}
}

bool match_order (const Match * l, const Match * r) {
	for (size_t i=0; i < l->size(); i++) {
		if (l->at(i) < r->at(i)) return true;
		if (l->at(i) > r->at(i)) return false;
	}
	return false;
}

void build_d2(ST_graph<Final_match> * d2_ptr, const vector<string> & strings, ST_graph<Match> * d1_ptr) {
	// Create D2 and add its source
	ST_graph<Match> & d1 = *d1_ptr;
	ST_graph<Final_match> & d2 = *d2_ptr;
	d2.source = d2.add_node(Final_match({&d1.get_source()}));

	queue<node_id_t> q;
	q.push(d2.source);

	for (; !q.empty(); q.pop()) {
		node_id_t curr_node = q.front();
		// Get the children and group them according to their character
		unordered_map<char, vector<Match*>> children; 
		for (const node_id_t m: d2.get_node(curr_node).modeled) {
			for (const node_id_t child_pos: d1.get_children_of(m)) {
				Match * child = &(d1.get_node(child_pos));
				children[strings[0][child->at(0)]].push_back(child); 
			}
		}
		// order of nodes is NOT deterministic and consistent
		// This is key for the unordered map to work with strings of Final_match

		vector<Match> virtual_matches;
		// build virtual matches: O(n*K) for each O(|sigma|) virt match + O(|sigma|^2)
		for (auto ch_iter = children.begin(); ch_iter != children.end(); ++ch_iter) {
			if (ch_iter->second.size() > 0) virtual_matches.push_back(Match (fishlet(ch_iter->second)));
		}

		// filter out matches that are dominated by some virtual match: O(n*sigma)*sigma?
		for (auto ch_iter = children.begin(); ch_iter != children.end(); ++ch_iter) {
			for (size_t i = 0; i < ch_iter->second.size(); i++) {
				Match & child = *(ch_iter->second[i]);
				for (size_t v = 0; v < virtual_matches.size(); ++v) {
					if (virtual_matches[v] < child) {
						ch_iter->second.erase(ch_iter->second.begin()+i); // TODO: do not use erase
						i--;
						break;
					}
				}
			}
		}
		
		// create d2 nodes with the matches that survived
		for (auto ch_iter = children.begin(); ch_iter != children.end(); ++ch_iter) {
			if (ch_iter->second.size() > 0) {
				// TODO O(n log n) is this a bottleneck??!!!
				sort(ch_iter->second.begin(), ch_iter->second.end(), match_order);
				Final_match final_child(ch_iter->second);
				if (!(d2.contains_match(final_child))) {
					if (DEBUG) cout <<"d2 does not contain "<<final_child<<endl;
					q.push(d2.add_node(final_child));
					if (ch_iter->second[0] == &(d1.get_sink()))
						d2.sink=d2.get_pos(final_child);
				}
				d2.add_edge(curr_node,d2.get_pos(final_child));
				if (DEBUG || d_flag) cout << d2.get_node(curr_node).to_string() <<" -> " << final_child.to_string() <<endl;
			}
		}
	}
			
}


// fill the rank&select data structure
void fill_preprocessing_occurrence(unordered_map<char,vector<vector<int>>> & occurrence, const vector<string> & strings, const vector<char> & sigma, bool previous) {
	int K = strings.size();
	for (char c: sigma) {
		occurrence[c] = vector<vector<int>>(K); 
		if (DEBUG) cout << c <<endl;
		for (int i = 0; i < K; i++) {
			if (DEBUG) cout << i <<endl;
			string curr_str = strings[i];
			occurrence[c][i] = vector<int>(curr_str.size(),-1);
			for (size_t j = 1; j < curr_str.size(); j++) {
				if (previous) {
					if (curr_str[j-1] == c)
						occurrence[c][i][j] = j-1;
					else
						occurrence[c][i][j] = occurrence[c][i][j-1];
					if (DEBUG) cout<<occurrence[c][i][j]<<" ";
				} else { // get the next occurrence by starting from the end
					int true_j = curr_str.size()-j-1;
					if (curr_str[true_j+1] == c)
						occurrence[c][i][true_j] = true_j+1;
					else
						occurrence[c][i][true_j] = occurrence[c][i][true_j+1];
					if (DEBUG) cout<<occurrence[c][i][true_j]<<" ";
				}
			}
			if (DEBUG) cout<<"\n";
		}
	}
}

int main(int argc, char* argv[]) {
	srand(time(NULL));
	int K = 2;
	bool m_flag = false;
	bool l_flag = false;
	bool a_flag = false;
	bool z_flag = false;
	bool r_flag = false;
	bool i_flag = false;
	int alphabet_len = 4;
	vector<string> strings;
	int len=15;

	for (int i = 1; i < argc; i++) {
		std::string arg = argv[i];
		if (arg == "-k" ) {
			if (i + 1 < argc)
				K = stoi(argv[i + 1]);
			else
				cout<< "flag -k must be followed by the requested number of strings\n";
			i++;
		} else if (arg == "-s" ) {
			if (i + 1 < argc)
				alphabet_len = stoi(argv[i + 1]);
			else
				cout<< "flag -s must be followed by the requested size of the alphabet\n";
			i++;
		} else if (arg == "-n" ) {
			if (i + 1 < argc)
				len = stoi(argv[i + 1]);
			else
				cout<< "flag -n must be followed by the requested length of the strings\n";
			i++;
		} else if (arg == "-m") {
			m_flag = true;
		} else if (arg == "-d") {
			d_flag = true;
		} else if (arg == "-i") {
			i_flag = true;
		} else if (arg == "-z") {
			z_flag = true;
		} else if (arg == "-a") {
			a_flag = true;
		} else if (arg == "-r") {
			r_flag = true;
		} else if (arg == "-l") {
			l_flag = true;
		} else if (arg == "-h" || arg == "--help") {
			cout << "usage:\n\t./mcdag [FLAG...] [STRING...]\n\n"
				<< "Flags:\n\t-k num_strings\tset the number of"
				<< "strings to compute mcs on\n"
				<< "\t-n len_strings\t\tset the length of the strings"
				<< "\n\t-s alphabet_size\tset the alphabet size\n"
				<< "\t-m\t\t\tprint all mcs (!!)\n"
				<< "\t-d\t\t\tprint McDag\n"
				<< "\t-l\t\t\tprint the distribution of the mcs lengths\n"
				<< "\t-r\t\t\tprint the number of lcs-1 and lcs\n"
				<< "\t-z\t\t\tminimalize McDag\n"
				<< "\t-a\t\t\tapply all flags to approximate indices\n"
				<< "\t-i\t\t\tdo not use McDag optimizations (slow)\n"
				<< "\t-h, --help\t\tshow this help\n";
			return 0;
		} else {
			strings.push_back(arg);
		}
	}
	if (strings.size() == 1) { cout << "please pass at least 2 strings (or none)\n"; return -1;}

	vector<char> sigma; 
	sigma.reserve(93);
	for (int i = 0; i <256; i++) {
		char c = static_cast<char>((i+65)%256);
		if (c != '$' && c != '#' && isprint(c))
			sigma.push_back(c);
	}
	assert(alphabet_len <= sigma.size());
	sigma.erase(sigma.begin()+alphabet_len,sigma.end());
	if (strings.size() == 0) {
		// randomly generate strings
		strings = vector<string>(K);
		for (int i = 0; i < K; i++) {
			strings[i].reserve(len+2);
			strings[i] += '#';
			for (int j = 0; j < len; j++) {
				// Warning: this is not uniform https://c-faq.com/lib/randrange.html
				strings[i] += sigma[rand() % (alphabet_len)];
			}
			strings[i] += '$';
			cout << strings[i]<<endl;
		}
	} else {
		// populate sigma according to passed strings
		K = strings.size();
		unordered_set<char> sigma_set;
		for (int i = 0; i < K; i++)
			for (size_t j = 0; j < strings[i].size(); j++)
				sigma_set.insert(strings[i][j]);
		if (sigma_set.find('#') != sigma_set.end() || sigma_set.find('$') != sigma_set.end()) {
			cout << "please do not use characters # and $\n";
			return -2;
		}
		sigma = vector<char>(sigma_set.begin(),sigma_set.end());
		for (int i = 0; i < K; i++) {
			strings[i].insert(0,"#");
			strings[i] += "$";
			cout << strings[i]<<endl;
		}
	}
	cout << "strings lengths: ";
	for (int i = 0; i < K; i++) {
		cout << strings[i].size()-2<<" ";
	}
	cout << endl;

	int num_matches = 2;
	for (size_t i = 0; i < sigma.size(); i++) {
		char c = sigma[i];
		int matches_of_c = 1;
		for (size_t j = 0; j < strings.size(); j++) {
			matches_of_c*= count(strings[j].cbegin(), strings[j].cend(), c);
		}
		num_matches += matches_of_c;
	}
	cout << "number of possible matches: "<< num_matches << endl;

	ST_graph<Match>* d1_ptr;
	ST_graph<Match>* d1_codet;
	ST_graph<Final_match> * d2 = new ST_graph<Final_match>(true,true, true);

	gettimeofday(&t_begin, 0);
	unordered_map<char,vector<vector<int>>> prev_occurrence;
	fill_preprocessing_occurrence(prev_occurrence, strings, sigma, true);

	if (!i_flag) {
		// preprocess helper data structure of size sigma x K x n
		unordered_map<char,vector<vector<int>>> next_occurrence;
		fill_preprocessing_occurrence(next_occurrence, strings, sigma, false);

		d1_ptr = new ST_graph<Match>(a_flag,true, true); 
		build_d1(d1_ptr,strings,next_occurrence,sigma, false);
		gettimeofday(&t_end, 0);
		cout << "time for d1r: " << time_elapsed(t_begin, t_end) << " seconds."<< endl;

		cout << "number of deterministic d1 nodes: "<< d1_ptr->get_num_nodes() << ", number of deterministic d1 edges: "<< d1_ptr->get_num_edges()<<endl;
		if (a_flag) {
			if (d1_ptr->print_stats(z_flag, m_flag, l_flag, r_flag, false, "deterministic d1", strings, sigma, true) == -1)
				return -1;
		}

		d1_codet = new ST_graph<Match>(true,a_flag, false);
		codeterminize_d1(d1_codet, d1_ptr,strings,prev_occurrence, sigma);
		gettimeofday(&t_end, 0);
		cout << "time for d1pp: " << time_elapsed(t_begin, t_end) << " seconds."<< endl;
		delete d1_ptr;
		cout << "number of codeterministic d1_intersection nodes: "<< d1_codet->get_num_nodes() << ", number of codeterministic d1_intersection edges: "<< d1_codet->get_num_edges()<<endl;
		if (a_flag) {
			if (d1_codet->print_stats(z_flag, m_flag, l_flag, r_flag, false, "codeterministic d1_intersection", strings, sigma, true) == -1)
				return -1;
		}

		build_d2(d2, strings, d1_codet);
		gettimeofday(&t_end, 0);
		cout << "time for d2 through d1pp : " << time_elapsed(t_begin, t_end) << " seconds."<< endl;
		delete d1_codet;
	} else {
		ST_graph<Match>* csa_ptr;
		csa_ptr = new ST_graph<Match>(true,a_flag, false); 
		build_csa(csa_ptr,strings,prev_occurrence,sigma, true);
		gettimeofday(&t_end, 0);
		cout << "time for csa: " << time_elapsed(t_begin, t_end) << " seconds."<< endl;
		cout << "number of deterministic csa nodes: "<< csa_ptr->get_num_nodes() << ", number of deterministic csa edges: "<< csa_ptr->get_num_edges()<<endl;
		if (a_flag) {
			if (csa_ptr->print_stats(z_flag, m_flag, l_flag, r_flag, false, "deterministic csa", strings, sigma, false) == -1)
				return -1;
		}
		// d1_ptr = new ST_graph<Match>(true,a_flag, false); 
		// build_d1(d1_ptr,strings,prev_occurrence,sigma);
		// cout << "number of deterministic d1 nodes: "<< d1_ptr->get_num_nodes() << ", number of deterministic d1 edges: "<< d1_ptr->get_num_edges()<<endl;
		// if (a_flag) {
		// 	if (d1_ptr->print_stats(z_flag, m_flag, l_flag, r_flag, false, "codeterministic d1", strings, sigma, true) == -1)
		// 		return -1;
		// }
		build_d2(d2, strings, csa_ptr);
		gettimeofday(&t_end, 0);
		cout << "time for d2 through csa : " << time_elapsed(t_begin, t_end) << " seconds."<< endl;
		// delete d1_ptr;
		delete csa_ptr;
	}
	
	cout << "number of d2 nodes: "<< d2->get_num_nodes()<< ", number of d2 edges: "<< d2->get_num_edges()<<endl;
	if (d2->print_stats(z_flag, m_flag, l_flag, r_flag, false, "d2", strings, sigma, true) == -1)
		return -1;

	delete d2;
	return 0;
}
