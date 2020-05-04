#include <stdio.h>
#include <stdlib.h>
#include "glpk.h"


void readGraph(int *N, int *M,
               int **from, int **to, int **weight) {
    scanf("%d %d", N, M);

    *from = (int *) calloc(*M, sizeof(int));
    *to = (int *) calloc(*M, sizeof(int));
    *weight = (int *) calloc(*M, sizeof(int));

    for (int i = 0; i < *M; ++i) {
        scanf("%d %d %d", (*from) + i, (*to) + i, (*weight) + i);
    }
}


void readGraphFile(int *N, int *M,
               int **from, int **to, int **weight,
               FILE *file) {
    fscanf(file, "%d %d", N, M);

    *from = (int *) calloc(*M, sizeof(int));
    *to = (int *) calloc(*M, sizeof(int));
    *weight = (int *) calloc(*M, sizeof(int));

    for (int i = 0; i < *M; ++i) {
        fscanf(file, "%d %d %d", (*from) + i, (*to) + i, (*weight) + i);
    }
}


void buildMatchingProb(glp_prob *problem, int N, int M,
                       int *from, int *to, int *weight) {

    glp_set_obj_dir(problem, GLP_MIN);

    // Add edges` variables
    int k;
    for (int i = 0; i < M; ++i) {
        k = glp_add_cols(problem, 1);

        glp_set_col_name(problem, k, "");
        glp_set_col_kind(problem, k, GLP_BV);
        glp_set_obj_coef(problem, k, weight[i]);
    }

    // Add vertex constraints
    int len;
    int *ind = (int *) calloc(1 + M, sizeof(int));
    double *val = (double *) calloc(1 + M, sizeof(double));
    for (int i = 0; i < N; ++i) {
        len = 0;
        for (int j = 0; j < M; ++j) {
            if ((from[j] != i) && (to[j] != i)) {
                continue;
            }

            ++len;
            ind[len] = j + 1;
            val[len] = 1.;
        }

        k = glp_add_rows(problem, 1);

        glp_set_row_name(problem, k, "");
        glp_set_row_bnds(problem, k, GLP_FX, 1., 1.);
        glp_set_mat_row(problem, k, len, ind, val);
    }

    free(ind);
    free(val);
}


void freeStoredData(int *from, int *to, int *weight) {
    free(from);
    free(to);
    free(weight);
}


int buildMatching(glp_prob *problem, int N, int M,
                  int *minMatchingWeight, int **edgesInd) {
    int status = glp_mip_status(problem);
    if ((status != GLP_FEAS) && (status != GLP_OPT)) {
        return 1;
    }

    *minMatchingWeight = glp_mip_obj_val(problem);

    *edgesInd = (int *) calloc(M, sizeof(int));
    for (int i = 0; i < M; ++i) {
        (*edgesInd)[i] = glp_mip_col_val(problem, 1 + i);
    }

    return 0;
}


void showMatching(int N, int M, int minMatchingWeight,
                  int *edgesInd) {
    printf("%d\n", minMatchingWeight);

    for (int i = 0; i < M; ++i) {
        if (edgesInd[i] == 0) {
            continue;
        }

        printf("%d ", i);
    }
    printf("\n");
}


int main(int argc, char *argv[]) {
    int N, M;
    int *from, *to;
    int *weight;

    if (argc > 1) {
        FILE *file;
        file = fopen(argv[1], "r");

        readGraphFile(&N, &M, &from, &to, &weight, file);
    } else {
        readGraph(&N, &M, &from, &to, &weight);
    }

    glp_prob *problem = glp_create_prob();
    glp_iocp iocp;

    buildMatchingProb(problem, N, M, from, to, weight);
    freeStoredData(from, to, weight);

    // MIP PARAMS INIT
    glp_init_iocp(&iocp);
    iocp.msg_lev = GLP_MSG_OFF;
    iocp.br_tech = GLP_BR_DTH; /* most fractional variable */
    iocp.bt_tech = GLP_BT_BPH; /* best local bound */
    iocp.sr_heur = GLP_OFF; /* disable simple rounding heuristic */
    iocp.gmi_cuts = GLP_ON; /* enable Gomory cuts */
    iocp.mir_cuts = GLP_ON;
    iocp.presolve = GLP_ON; /* LP-relaxation pre-solution */

    glp_intopt(problem, &iocp);

    int minMatchingWeight;
    int *edgesInd;
    if (buildMatching(problem, N, M, &minMatchingWeight, &edgesInd)) {
        glp_delete_prob(problem);
        printf("Optimal solution undefined\n");
        return 0;
    }
    showMatching(N, M, minMatchingWeight, edgesInd);

    free(edgesInd);
    glp_delete_prob(problem);
}
