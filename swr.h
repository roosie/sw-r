//Roosie implemention of Smith-Waterman 2020
//TODO implement affine gap penalty function

#include <iostream>
#include <fstream>
#include <string>
#include <bits/stdc++.h>

using namespace std;

struct DP_cell{
    int score;
    int gap_count;
    int parent_i;
    int parent_j;
    bool visited = false;
};

struct point{
    int i;
    int j;
};

class swr{
    public:
        swr(string new_a_sequence, string new_b_sequence, int new_match, int new_mismatch);
        int score_calc(int i, int j);
        int gap_penalty(int i, int j);
        int compute_grid();
        int print_grid();
        int hello_world();
        int max(int a, int b, int c, int d);
        string generate_spaces_for_output(int number, int spaces);
        point largest_alignment_index(vector<point*> set);
        vector<point*> determine_largest();
        void print_local_alignment(point meeting);
    private: 
        vector<vector<DP_cell*>> grid;
        string a_sequence;
        string b_sequence;
        int match;
        int mismatch;

};