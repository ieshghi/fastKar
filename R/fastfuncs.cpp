// Fast versions of the graph functions 
#include <Rcpp.h>
#include <vector>
#include <set>
#include <algorithm>
using namespace Rcpp;

// [[Rcpp::export]]
List traverse_graph_cpp(DataFrame A, NumericVector loose_ends) {
  // Extract columns from the data frame
  IntegerVector n = A["n"];
  IntegerVector left = A["left"];
  IntegerVector right = A["right"];
  IntegerVector id = A["id"];
  
  int nrow = n.size();
  
  // Data structures for tracking
  std::set<int> visited_edges;
  std::set<int> visited_rows;
  List traversed_paths;
  
  // Make a copy of loose_ends to modify
  std::set<int> remaining_loose_ends;
  for (int i = 0; i < loose_ends.size(); i++) {
    remaining_loose_ends.insert(loose_ends[i]);
  }
  
  // Remove loose ends that have already been visited
  for (auto it = remaining_loose_ends.begin(); it != remaining_loose_ends.end();) {
    if (visited_edges.find(*it) != visited_edges.end()) {
      it = remaining_loose_ends.erase(it);
    } else {
      ++it;
    }
  }
  
  // First part: traverse paths starting from loose ends
  while (!remaining_loose_ends.empty()) {
    int start_edge = *remaining_loose_ends.begin();  // Pick an unvisited loose end
    std::vector<int> path;
    std::vector<int> nodepath;
    int current_edge = start_edge;
    
    while (true) {
      path.push_back(current_edge);
      visited_edges.insert(current_edge);
      
      // Find the row where current_edge appears in 'left' or 'right'
      int row_idx = -1;
      for (int i = 0; i < nrow; i++) {
        if ((left[i] == current_edge || right[i] == current_edge) && 
            visited_rows.find(id[i]) == visited_rows.end()) {
          row_idx = i;
          break;
        }
      }
      
      if (row_idx == -1) break;  // No valid row found
      
      visited_rows.insert(id[row_idx]);
      
      // Determine the next edge
      int next_edge;
      if (left[row_idx] == current_edge) {
        nodepath.push_back(n[row_idx]);
        next_edge = right[row_idx];
      } else {
        nodepath.push_back(-n[row_idx]);
        next_edge = left[row_idx];
      }
      
      // Check if next_edge is a loose end and we've already traversed something
      if (!path.empty() && 
          remaining_loose_ends.find(next_edge) != remaining_loose_ends.end()) {
        break;
      }
      
      // Find the next row containing next_edge
      int next_row_idx = -1;
      for (int i = 0; i < nrow; i++) {
        if ((left[i] == next_edge || right[i] == next_edge) && 
            visited_rows.find(id[i]) == visited_rows.end()) {
          next_row_idx = i;
          break;
        }
      }
      
      if (next_row_idx == -1) {
        // No more paths, add the final edge
        visited_edges.insert(next_edge);
        path.push_back(next_edge);
        break;
      }
      
      current_edge = next_edge;  // Move to the next edge
    }
    
    // Add the completed path
    traversed_paths.push_back(nodepath);
    
    // Update unvisited loose ends
    remaining_loose_ends.clear();
    for (int i = 0; i < loose_ends.size(); i++) {
      if (visited_edges.find(loose_ends[i]) == visited_edges.end()) {
        remaining_loose_ends.insert(loose_ends[i]);
      }
    }
  }
  
  // Second part: detect cycles
  List cycles;
  
  // Find remaining unvisited rows
  std::vector<int> remaining_row_indices;
  for (int i = 0; i < nrow; i++) {
    if (visited_rows.find(id[i]) == visited_rows.end()) {
      remaining_row_indices.push_back(i);
    }
  }
  
  while (!remaining_row_indices.empty()) {
    int start_idx = remaining_row_indices[0];
    int start_edge = right[start_idx];  // Start from right edge
    
    visited_rows.insert(id[start_idx]);
    
    std::vector<int> path;
    std::vector<int> nodepath;
    nodepath.push_back(n[start_idx]);  // Start with the node value
    
    int current_edge = start_edge;
    
    while (true) {
      path.push_back(current_edge);
      
      // Find the row containing this edge
      int row_idx = -1;
      for (int i : remaining_row_indices) {
        if ((left[i] == current_edge || right[i] == current_edge) && 
            visited_rows.find(id[i]) == visited_rows.end()) {
          row_idx = i;
          break;
        }
      }
      
      if (row_idx == -1) {
        // No more valid rows, cycle is complete
        cycles.push_back(nodepath);
        break;
      }
      
      visited_rows.insert(id[row_idx]);
      
      // Determine the next edge
      int next_edge;
      if (left[row_idx] == current_edge) {
        nodepath.push_back(n[row_idx]);
        next_edge = right[row_idx];
      } else {
        nodepath.push_back(-n[row_idx]);
        next_edge = left[row_idx];
      }
      
      current_edge = next_edge;  // Move to next edge
    }
    
    // Update remaining row indices
    remaining_row_indices.clear();
    for (int i = 0; i < nrow; i++) {
      if (visited_rows.find(id[i]) == visited_rows.end()) {
        remaining_row_indices.push_back(i);
      }
    }
  }
  
  // Construct and return the result list
  List result;
  result["paths"] = traversed_paths;
  result["cycles"] = cycles;
  return result;
}