#include "swr.h"

int swr::hello_world(){
    cout << "hellow wolrd" << endl;
    return 0;
}

swr::swr(string new_a_sequence, string new_b_sequence, int new_match, int new_mismatch){
    cout << "swr constructor" << endl;
    //create a grid of proper size
    vector<vector<DP_cell*>> temp_grid(new_a_sequence.size()+1, vector<DP_cell*>(new_b_sequence.size()+1));
    //set the first row to 0
    
    for(int i = 0; i < temp_grid.size(); i++){
        temp_grid[i][0] = new DP_cell();
        temp_grid[i][0]->score = 0;
    }    //set the first column to 0
    for(int j = 0; j < temp_grid[0].size(); j++){
        temp_grid[0][j] = new DP_cell();
        temp_grid[0][j]->score = 0;
    }

    grid = temp_grid;
    match = new_match;
    mismatch = new_mismatch;
    a_sequence = new_a_sequence;
    b_sequence = new_b_sequence;
    
}

int swr::score_calc(int i, int j){
    //defaults to a match for upper left value
    int upper_left = grid[i-1][j-1]->score + match;
    //need to refactor gap_penalty to use affine gap penalty function 
    int left = grid[i-1][j]->score + gap_penalty(i,j);
    int upper = grid[i][j-1]->score + gap_penalty(i,j);

    //sets upper_left to a mismatch value if a mismatch occurs 
    if (a_sequence[i-1] != b_sequence[j-1]){
        upper_left -= match;
        upper_left += mismatch;
    }

    //creates cell
    grid[i][j] = new DP_cell();
    //cout << match << " " <<  mismatch << endl;
    grid[i][j]->score = max(upper, upper_left, left, 0);
    if(grid[i][j]->score == upper){
        grid[i][j]->parent_i = i;
        grid[i][j]->parent_j = j-1;
    }
    else if(grid[i][j]->score == left){
        grid[i][j]->parent_i = i-1;
        grid[i][j]->parent_j = j;
    }
    else {
        grid[i][j]->parent_i = i-1;
        grid[i][j]->parent_j = j-1;
    }
    return 1;
}

//should update to use afnity gap penalty
int swr::gap_penalty(int i, int j){
    return -5;
}

//finds max of 4 ints 
int swr::max(int a, int b, int c, int d){
    int e = a;
    if(a < b){
        e = b;
    }
    int f = c;
    if(c < d){
        f = d;
    }
    if(f > e){
        return f;
    }
    return e;
}

int swr::compute_grid(){
    cout << "compute grid" << endl;
    //calls calculate score on each node 
    
    for (int j = 1; j < grid[0].size(); j++){   
        for(int i = 1; i < grid.size(); i++){
            score_calc(i, j);
        }
    }
    return 1;
}

//prints grid
int swr::print_grid(){
    cout << "print grid" << endl;
    for(int j = 0; j < grid[0].size(); j++){
        for(int i = 0; i < grid.size(); i++){
        
            cout << grid[i][j]->score << " ";
        }
        cout << endl;
    }
    return 1;
}

//linear select to group the largest values 
vector<point*> swr::determine_largest(){
    cout << "determine largest" << endl;
    vector<point*> current_max;
    point* current;
    int max_value = 0;
    for(int j = 0; j < grid[0].size(); j++){
        for(int i = 0; i < grid.size(); i++){
            if(grid[i][j]->score > max_value){
                current_max.clear();
                current = new point();
                current->i = i;
                current->j = j;
                current_max.push_back(current);
                max_value = grid[i][j]->score;
            }
            else if(grid[i][j]->score == max_value){
                current = new point();
                current->i = i;
                current->j = j;
                current_max.push_back(current);
            }
        }
    }
    return current_max;
}

//consider just keep track of length during initial grid construction  
point swr::largest_alignment_index(vector<point*> set){
    cout << "largest algiment index" << endl;
    int max_len = 0;
    int current_length = 0;
    point max_set;
    int i = 0; 
    int j = 0;
    DP_cell* node;
    for(int index = 0; index < set.size(); index++){
        i = set[set.size() - index - 1]->i;
        j = set[set.size() - index - 1]->j;
        current_length = 0;
        node = grid[i][j];
        //tranverses parents node while counting distance travel, sets node to visited so theyre not visited again
        while(node->score > 0 && node->visited == false){
            i = node->parent_i;
            j = node->parent_j;
            current_length += 1;
            node->visited = true;
            node = grid[i][j];
        }
        //sets new max
        if(current_length > max_len){
            max_set.i = set[set.size() - index - 1]->i;
            max_set.j = set[set.size() - index - 1]->j;
            max_len = current_length;
        }
    }

    return max_set;
}
//assuminge 80x24 terminal window 
void swr::print_local_alignment(point meeting){
    //prints initial stats
    int a_len = a_sequence.size();
    int b_len = b_sequence.size();
    cout << "Scores: match " << match << ", mismatch " << mismatch << ", indel -5" << endl << endl;
    cout << "Sequence 1 = \"s1\", length = " << a_len << " characters" << endl; 
    cout << "Sequence 2 = \"s2\", length = " << b_len << " characters" << endl; 

    //calculates the offset distance
    int diff = meeting.i - meeting.j;
    int a_offset = 0;
    int b_offset = 0;
    if(diff < 0){
        a_offset += diff;
    }
    else{
        b_offset -= diff;
    }

    int max_dist = max(a_len - a_offset,b_len - b_offset, 0, 0);
    int index = 0;
    int temp_index = 0;
    int matches = 0;
    int gap_count = 0;
    int gap_running = 0;
    cout << a_offset << "," << b_offset << endl;

    //walks strings, counts number of gaps and number of matches 
    while(index < max_dist){
        cout << "s1 " << generate_spaces_for_output(index+1,4) << index+1 << "    ";
        temp_index = index;
        while(temp_index < index+60){
            if(temp_index+a_offset < 0 || temp_index+a_offset >= a_len){
                cout << "-";
            }
            else{
                cout << a_sequence[temp_index+a_offset];
            }
            temp_index += 1;
        }
        cout << endl;
        cout << "   " << "    " << "    ";
        temp_index = index;
        while(temp_index < index+60){
            if((temp_index+a_offset < 0 || temp_index+a_offset >= a_len) || (temp_index+b_offset < 0 || temp_index+b_offset >= b_len)){
                cout << " ";
                gap_running += 1;
                if(gap_running == 1){
                    gap_count +=1;
                }
            }
            else if (a_sequence[temp_index+a_offset] != b_sequence[temp_index+b_offset]){
                cout << " ";
                gap_running += 1;
                if(gap_running == 1){
                    gap_count +=1;
                }
            }
            else{
                cout << "|";
                matches +=1;
                gap_running = 0;
            }
            temp_index+=1;
        }
        cout << endl;
        cout << "s2 " << generate_spaces_for_output(index+1,4) << index+1 << "    ";
        temp_index = index;
        while(temp_index < index+60){
            if(temp_index+b_offset < 0 || temp_index+b_offset >= b_len){
                cout << "-";
            }
            else{
                cout << b_sequence[temp_index+b_offset];
            }
            temp_index += 1;
        }
        cout << endl;
        index += 60;
    }
    //prints result data
    cout << "Number of: matches = " << matches << ", mismatches = " << max_dist - matches << ", gaps = " << gap_count << endl;

}

string swr::generate_spaces_for_output(int number, int spaces){
    int count = 0;
    int temp = 0;
    temp = number;
    while(temp != 0) {
        count++;
        temp /= 10;
    }
    string base = "";
    for(int i = 0; i < spaces-count; i++){
        base += ' ';
    }
    return base;
}