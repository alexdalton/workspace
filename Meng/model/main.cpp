//  Copyright 2013 Google Inc. All Rights Reserved.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <gsl/gsl_rng.h>
#include "linelib.h"
#include "ransampl.h"

#define MAX_PATH_LENGTH 100

char node_file[MAX_STRING], link_file[MAX_STRING], output_file[MAX_STRING];
int binary = 0, num_threads = 1, vector_size = 100, negative = 5, iters = 10, epoch, mode = 0, model = 0, depth = 2;
long long samples = 1, edge_count_actual;
real alpha = 0.025, starting_alpha, restart = 0.3;

const gsl_rng_type * gsl_T;
gsl_rng * gsl_r;

line_node node0, node1;
line_hin hin;
line_trainer_edge trainer_edge1, trainer_edge2, trainer_edge3, trainer_edge4, trainer_edge5;
line_trainer_path trainer_path;

double func_rand_num()
{
    return gsl_rng_uniform(gsl_r);
}

void *training_thread(void *id)
{
    long long edge_count = 0, last_edge_count = 0;
    unsigned long long next_random = (long long)id;
    real *error_vec = (real *)calloc(vector_size, sizeof(real));
    real *error_p = (real *)calloc(vector_size, sizeof(real));
    real *error_q = (real *)calloc(vector_size, sizeof(real));
    int *node_lst = (int *)malloc(MAX_PATH_LENGTH * sizeof(int));

    while (1)
    {
        //judge for exit
        if (edge_count > samples / num_threads + 2) break;

        if (edge_count - last_edge_count > 1000)
        {
            edge_count_actual += edge_count - last_edge_count;
            last_edge_count = edge_count;
            printf("%cEpoch: %d/%d Alpha: %f Progress: %.3lf%%", 13, epoch + 1, iters, alpha, (real)edge_count_actual / (real)(samples + 1) * 100);
            fflush(stdout);
            alpha = starting_alpha * (1 - edge_count_actual / (real)(samples * iters+ 1));
            if (alpha < starting_alpha * 0.0001) alpha = starting_alpha * 0.0001;
        }

        trainer_edge1.train_sample(mode, alpha, error_vec, error_p, error_q, func_rand_num, next_random);
        trainer_edge2.train_sample(mode, alpha, error_vec, error_p, error_q, func_rand_num, next_random);
        trainer_edge3.train_sample(mode, alpha, error_vec, error_p, error_q, func_rand_num, next_random);
        trainer_edge4.train_sample(mode, alpha, error_vec, error_p, error_q, func_rand_num, next_random);
        trainer_edge5.train_sample(mode, alpha, error_vec, error_p, error_q, func_rand_num, next_random);

        edge_count += 2;
    }
    free(node_lst);
    free(error_vec);
    free(error_p);
    free(error_q);
    pthread_exit(NULL);
}

void TrainModel() {
    long a;
    pthread_t *pt = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
    starting_alpha = alpha;

    gsl_rng_env_setup();
    gsl_T = gsl_rng_rand48;
    gsl_r = gsl_rng_alloc(gsl_T);
    gsl_rng_set(gsl_r, 314159265);

    node0.init(node_file, vector_size);
    node1.init(node_file, vector_size);
    hin.init(link_file, &node0, &node1);
    trainer_edge1.init('a', &hin, negative);
    trainer_edge2.init('b', &hin, negative);
    trainer_edge3.init('c', &hin, negative);
    trainer_edge4.init('d', &hin, negative);
    trainer_edge5.init('e', &hin, negative);


    clock_t start = clock();
    printf("Training:");
    for (epoch = 0; epoch != iters; epoch++)
    {
        edge_count_actual = samples * epoch;
        mode = 0;
        for (a = 0; a < num_threads; a++) pthread_create(&pt[a], NULL, training_thread, (void *)a);
        for (a = 0; a < num_threads; a++) pthread_join(pt[a], NULL);

        edge_count_actual = samples * epoch;
        mode = 1;
        for (a = 0; a < num_threads; a++) pthread_create(&pt[a], NULL, training_thread, (void *)a);
        for (a = 0; a < num_threads; a++) pthread_join(pt[a], NULL);
    }
    printf("\n");
    clock_t finish = clock();
    printf("Total time: %lf\n", (double)(finish - start) / CLOCKS_PER_SEC);

    node0.output(output_file, binary);
    node1.output(output_file, binary);
}

int ArgPos(char *str, int argc, char **argv) {
    int a;
    for (a = 1; a < argc; a++) if (!strcmp(str, argv[a])) {
        if (a == argc - 1) {
            printf("Argument missing for %s\n", str);
            exit(1);
        }
        return a;
    }
    return -1;
}

int main(int argc, char **argv) {
    int i;
    if (argc == 1) {
        printf("HIN2VEC\n\n");
        printf("Options:\n");
        printf("Parameters for training:\n");
        printf("\t-node <file>\n");
        printf("\t\tA dictionary of all nodes\n");
        printf("\t-link <file>\n");
        printf("\t\tAll links between nodes. Links are directed.\n");
        printf("\t-path <int>\n");
        printf("\t\tAll meta-paths. One path per line.\n");
        printf("\t-output <int>\n");
        printf("\t\tThe output file.\n");
        printf("\t-binary <int>\n");
        printf("\t\tSave the resulting vectors in binary moded; default is 0 (off)\n");
        printf("\t-size <int>\n");
        printf("\t\tSet size of word vectors; default is 100\n");
        printf("\t-negative <int>\n");
        printf("\t\tNumber of negative examples; default is 5, common values are 5 - 10 (0 = not used)\n");
        printf("\t-samples <int>\n");
        printf("\t\tSet the number of training samples as <int>Million\n");
        printf("\t-iters <int>\n");
        printf("\t\tSet the number of interations.\n");
        printf("\t-threads <int>\n");
        printf("\t\tUse <int> threads (default 1)\n");
        printf("\t-alpha <float>\n");
        printf("\t\tSet the starting learning rate; default is 0.025\n");
        printf("\nExamples:\n");
        printf("./hin2vec -node node.txt -link link.txt -path path.txt -output vec.emb -binary 1 -size 100 -negative 5 -samples 5 -iters 20 -threads 12\n\n");
        return 0;
    }
    output_file[0] = 0;
    if ((i = ArgPos((char *)"-node", argc, argv)) > 0) strcpy(node_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-link", argc, argv)) > 0) strcpy(link_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-output", argc, argv)) > 0) strcpy(output_file, argv[i + 1]);
    if ((i = ArgPos((char *)"-binary", argc, argv)) > 0) binary = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-size", argc, argv)) > 0) vector_size = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-model", argc, argv)) > 0) model = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-negative", argc, argv)) > 0) negative = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-samples", argc, argv)) > 0) samples = (long long)(atof(argv[i + 1])*1000000);
    if ((i = ArgPos((char *)"-iters", argc, argv)) > 0) iters = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-depth", argc, argv)) > 0) depth = atoi(argv[i + 1]);
    if ((i = ArgPos((char *)"-restart", argc, argv)) > 0) restart = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-alpha", argc, argv)) > 0) alpha = atof(argv[i + 1]);
    if ((i = ArgPos((char *)"-threads", argc, argv)) > 0) num_threads = atoi(argv[i + 1]);
    TrainModel();
    return 0;
}