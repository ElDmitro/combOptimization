#include "glpk.h"

int main() {
    glp_prob *lp = glp_create_prob();
    glp_delete_prob(lp);
}
