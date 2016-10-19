#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <Eigen/Dense>
#include "ransampl.h"
#include <iostream>

#define MAX_STRING 500
#define EXP_TABLE_SIZE 1000
#define MAX_EXP 6
const int neg_table_size = 1e8;
const int hash_table_size = 30000000;

typedef float real;

typedef Eigen::Matrix< real, Eigen::Dynamic,
Eigen::Dynamic, Eigen::RowMajor | Eigen::AutoAlign >
BLPMatrix;

typedef Eigen::Matrix< real, 1, Eigen::Dynamic,
Eigen::RowMajor | Eigen::AutoAlign >
BLPVector;

struct struct_node {
    char *word;
    char type;
};

struct hin_nb {
    int nb_id;
    double eg_wei;
    char eg_tp;
};

class line_node
{
protected:
    struct struct_node *node;
    int node_size, node_max_size, vector_size;
    char node_file[MAX_STRING];
    int *node_hash;
    real *_vec;
    Eigen::Map<BLPMatrix> vec;
    
    int get_hash(char *word);
    int add_node(char *word, char type);
public:
    line_node();
    ~line_node();
    
    friend class line_link;
    friend class line_hin;
    friend class line_trainer_edge;
    friend class line_trainer_path;
    
    void init(char *file_name, int vector_dim);
    int search(char *word);
    void output(char *file_name, int binary);
    
    //friend void linelib_output_batch(char *file_name, int binary, line_node **array_line_node, int cnt);
};

class line_hin
{
protected:
    char hin_file[MAX_STRING];
    
    line_node *node_u, *node_v;
    std::vector<hin_nb> *hin;
    long long hin_size;
    
public:
    line_hin();
    ~line_hin();
    
    friend class line_trainer_edge;
    friend class line_trainer_path;
    
    void init(char *file_name, line_node *p_u, line_node *p_v);
};

class line_trans
{
protected:
    std::string type;
    int vector_size;
    real *_P, *_Q;
    real bias;
    Eigen::Map<BLPVector> P, Q;
public:
    line_trans();
    ~line_trans();
    
    friend class line_trainer_edge;
    friend class line_trainer_path;
    
    void init(std::string map_type, int vector_dim);
};

class line_trainer_edge
{
protected:
    line_hin *phin;
    
    int *u_nb_cnt; int **u_nb_id; double **u_nb_wei;
    double *u_wei, *v_wei;
    ransampl_ws *smp_u, **smp_u_nb;
    real *expTable;
    int neg_samples, *neg_table;
    
    char edge_tp;
    
    static std::vector<line_trans *> vec_trans;
    static std::map<std::string, int> map_trans;
    static int cnt_trans;
public:
    line_trainer_edge();
    ~line_trainer_edge();
    
    static void print_trans();
    
    void init(char edge_type, line_hin *p_hin, int negative);
    void train_sample(int mode, real alpha, real *_error_vec, real *_error_p, real *_error_q, double (*func_rand_num)(), unsigned long long &rand_index);
    void train_sample_depth(int mode, real alpha, int depth, real *_error_vec, real *_error_p, real *_error_q, double (*func_rand_num)(), unsigned long long &rand_index);
    void train_sample_randwalk(int mode, real alpha, real restart, real *_error_vec, real *_error_p, real *_error_q, double (*func_rand_num)(), unsigned long long &rand_index);
};

class line_trainer_path
{
protected:
    line_hin *phin;
    
    real *expTable;
    int neg_samples;
    double **dp_cnt;
    int **neg_table;
    ransampl_ws *smp_init;
    double *smp_init_weight;
    int *smp_init_index;
    ransampl_ws ***smp_dp;
    double ***smp_dp_weight;
    int ***smp_dp_index;
    
    std::string path, path_node, path_link;
    int path_size;
    
    static std::vector<line_trans *> vec_trans;
    static std::map<std::string, int> map_trans;
    static int cnt_trans;
    
    void sample_path(int *node_lst, double (*func_rand_num)());
public:
    line_trainer_path();
    ~line_trainer_path();
    
    static void print_trans();
    
    int get_path_length();
    void init(std::string meta_path, line_hin *p_hin, int negative);
    void train_path(int mode, int *node_lst, real alpha, real *_error_vec, real *_error_p, real *_error_q, double (*func_rand_num)(), unsigned long long &rand_index, int model);
};

