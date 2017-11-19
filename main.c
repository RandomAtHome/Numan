#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#define _MAIN_ELEM

double DETERMINANT = 1;
/* operation structure block*/
enum OPERATIONS {
    NONE = 0,
    SWAP = 1,
    SUBTRACT = 2,
    SUBTRACT_B = 3,
    NORMALIZE = 4
};

typedef struct operation {
    int type;
    size_t which;
    size_t from;
    double mod;
} Operation;

void add_operation(Operation **ops, size_t *op_len, size_t *op_i, int type, size_t from, size_t which, double mod) {
    if (!ops) {
        return;
    }
    if (*op_i == *op_len) {
        *op_len *= 2;
        *ops = realloc(*ops, *op_len * sizeof(Operation));
    }
    (*ops)[*op_i].type = type;
    (*ops)[*op_i].which = which;
    (*ops)[*op_i].from = from;
    (*ops)[*op_i].mod = mod;
    (*op_i)++;
    return;
}

/*input and matrix creation block*/
void input_fill_matrix(double **matrix, size_t side) {
    size_t i, j;
    printf("Enter matrix:\n");
    for (i = 0; i < side; i++) {
        for (j = 0; j < side; j++) {
            scanf("%lf", &matrix[i][j]);
        }
    }
    printf("Enter values vector:\n");
    for (i = 0; i < side; i++) {
        scanf("%lf", &matrix[i][side]);
    }
    return;
}

void formula_filled_matrix(double **matrix, size_t side) {
    printf("Enter x for formula\n");
    double x;
    scanf("%lf", &x);
    double q = 0.989;
    size_t i, j;
    for (i = 0; i < side; i++) {
        for (j = 0; j < side; j++) {
            if (j != i) {
                matrix[i][j] = pow(q, i + j + 2) + 0.1 * ((int) (j - i));
            } else {
                matrix[i][j] = pow(q - 1, i + j + 2);
            }
        }
    }
    for (i = 0; i < side; i++) {
        matrix[i][side] = x * exp(x / (i + 1)) * cos(x / (i + 1));
    }
    return;
}

double **get_matrix(size_t *side) {
    printf("Use formulated matrix? (Y/n)\n");
    char *temp = NULL;
    char flag = 0;
    size_t len = 1;
    while (1) {
        getline(&temp, &len, stdin);
        if (temp[0] == '\n' || temp[0] == 'y') {
            flag = 1;
            break;
        } else if (temp[0] == 'n') {
            break;
        }
    }
    free(temp);
    if (!flag) {
        printf("Enter n - side of matrix:\n");
        scanf("%zu", side);
    } else {
        *side = 100;
    }
    double **matrix = calloc(*side, sizeof(double *));
    matrix[0] = calloc(*side * (*side + 1), sizeof(double));
    size_t i;
    for (i = 1; i < *side; i++) {
        matrix[i] = matrix[0] + i * (*side + 1);
    }
    if (!flag) {
        input_fill_matrix(matrix, *side);
    } else {
        formula_filled_matrix(matrix, *side);
    }
    return matrix;
}

double **get_identity_matrix(size_t side) {
    double **matrix = calloc(side, sizeof(double *));
    matrix[0] = calloc(side * (side + 1), sizeof(double));
    size_t i;
    for (i = 0; i < side; i++) {
        matrix[i] = matrix[0] + i * (side + 1);
        matrix[i][i] = 1;
    }
    return matrix;
}

/*printing block*/
void __print_matrix(double **matrix, size_t side) {
    size_t i, j;
    char var_name[13] = {0};
    printf("|");
    for (j = 0; j < side; j++) {
        snprintf(var_name, 13, "%s%zu", "x", j + 1);
        printf("%12s|", var_name);
    }
    printf("%12s|\n", "b");
    for (i = 0; i < side; i++) {
        printf("|");
        for (j = 0; j < side + 1; ++j) {
            printf("%12lg|", matrix[i][j]);
        }
        printf("\n");
    }
    return;
}

void __print_answers(double **matrix, size_t side) {
    size_t i;
    char var_name[13] = {0};
    for (i = 0; i < side; i++) {
        snprintf(var_name, 13, "%s%zu", "x", i + 1);
        printf("|%6s|%12lg|\n", var_name, matrix[i][side]);
    }
    return;
}

/*line normalization op block*/
int normalize_line(double **matrix, size_t side, size_t which, double *mod) {
    size_t j, first_elem = (size_t) -1;
    for (j = 0; j < side; j++) {
        if (matrix[which][j]) {
            first_elem = j;
            break;
        }
    }
    if (first_elem == (size_t) -1) {
        return -1;
    }
    *mod = matrix[which][first_elem];
    for (j = first_elem + 1; j < side + 1; j++) {
        matrix[which][j] /= *mod;
    }
    DETERMINANT *= matrix[which][first_elem];
    matrix[which][first_elem] = 1;
    return 0;
}

int r_normalize_line(double **matrix, size_t side, size_t which, double mod) {
    size_t j;
    for (j = 0; j < side + 1; j++) {
        matrix[which][j] /= mod;
    }
    return 0;
}

/*swap op block*/
int swap_lines(double **matrix, size_t side, size_t line1, size_t line2) {
    if (line1 == line2) {
        return 1;
    }
    double *temp_buf = malloc((side + 1) * sizeof(double));
    memcpy(temp_buf, matrix[line1], (side + 1) * sizeof(double));
    memcpy(matrix[line1], matrix[line2], (side + 1) * sizeof(double));
    memcpy(matrix[line2], temp_buf, (side + 1) * sizeof(double));
    free(temp_buf);
    DETERMINANT *= -1;
    return 0;
}

/*subtract ops block*/
size_t find_fitting_row(double **matrix, size_t side, size_t column) {
    size_t j, res = (matrix[column][column] == 0 ? (size_t) -1 : column);
#ifdef _MAIN_ELEM
    double max_val = fabs(matrix[column][column]);
#endif
    for (j = column; j < side; j++) {
#ifdef _MAIN_ELEM
        if (fabs(matrix[j][column]) > max_val) {
            res = j;
        }
#else
        if (matrix[j][column]) {
            return j;
        }
#endif
    }
    return res;
}

int r_subtract_line_from_line(double **matrix, size_t side, size_t from, size_t which, double mod) {
    size_t j, first_elem = (size_t) -1;
    for (j = 0; j < side; j++) {
        if (matrix[which][j]) {
            first_elem = j;
            break;
        }
    }
    if (first_elem == (size_t) -1) {
        return -1;
    }
    if (which != from) {
        for (j = first_elem; j < side + 1; j++) {
            matrix[from][j] -= mod * matrix[which][j];
        }
    }
    return 0;
}

int subtract_line_from_line(double **matrix, size_t side, size_t from, size_t which, double *mod) {
    size_t j, first_elem = (size_t) -1;
    for (j = 0; j < side; j++) {
        if (matrix[which][j]) {
            first_elem = j;
            break;
        }
    }
    if (first_elem == (size_t) -1) {
        return -1;
    }
    if (which != from) {
        *mod = matrix[from][first_elem] / matrix[which][first_elem];
        for (j = first_elem + 1; j < side + 1; j++) {
            matrix[from][j] -= *mod * matrix[which][j];
        }
        matrix[from][first_elem] = 0;
    }
    return 0;
}

int subtract_line_from_line_b(double **matrix, size_t side, size_t from, size_t which, double *mod) {
    if (!matrix[which][which]) {
        return 1;
    }
    *mod = matrix[from][which];
    matrix[from][side] -= matrix[which][side] * matrix[from][which];
    matrix[from][which] = 0;
    return 0;
}

/*eliminations functions block*/
int straight_move(double **matrix, size_t side, Operation **ops, size_t *op_len, size_t *op_i) {
    size_t column, row, s_row;
    double mod;
    for (column = 0; column < side; column++) {
        if ((s_row = find_fitting_row(matrix, side, column)) == (size_t) -1) {
            return 1;
        }
        for (row = column; row < side; row++) {
            if ((subtract_line_from_line(matrix, side, row, s_row, &mod) == -1)) {
                return 1;
            }
            add_operation(ops, op_len, op_i, SUBTRACT, row, s_row, mod);
        }
        swap_lines(matrix, side, s_row, column);
        add_operation(ops, op_len, op_i, SWAP, s_row, column, 0);
        normalize_line(matrix, side, column, &mod);
        add_operation(ops, op_len, op_i, NORMALIZE, (size_t) -1, column, mod);
    }
    return 0;
}

int backwards_move(double **matrix, size_t side, Operation **ops, size_t *op_len, size_t *op_i) {
    size_t s_row, row;
    double mod;
    for (s_row = side - 1; s_row; s_row--) {
        for (row = s_row - 1; row != (size_t) -1; row--) {
            if (subtract_line_from_line_b(matrix, side, row, s_row, &mod)) {
                return 1;
            }
            add_operation(ops, op_len, op_i, SUBTRACT_B, row, s_row, mod);
        }
    }
    return 0;
}

/*inverse matrix related func*/
int apply_operations(double **matrix, size_t side, Operation *ops, size_t op_len) {
    size_t op_i;
    for (op_i = 0; op_i < op_len; op_i++) {
        switch (ops[op_i].type) {
            case NONE:
                break;
            case SWAP:
                swap_lines(matrix, side, ops[op_i].from, ops[op_i].which);
                break;
            case SUBTRACT:
                if (r_subtract_line_from_line(matrix, side, ops[op_i].from, ops[op_i].which, ops[op_i].mod)) {
                    return SUBTRACT;
                }
                break;
            case SUBTRACT_B:
                if (r_subtract_line_from_line(matrix, side, ops[op_i].from, ops[op_i].which, ops[op_i].mod)) {
                    return SUBTRACT_B;
                }
                break;
            case NORMALIZE:
                r_normalize_line(matrix, side, ops[op_i].which, ops[op_i].mod);
            default:
                break;
        }
    }
    return 0;
}

/*main*/
int main() {
#ifdef _MAIN_ELEM
    printf("Compiled with _MAIN_ELEM\n");
#else
    printf("Compiled without _MAIN_ELEM\n");
#endif
    size_t side;
    double **matrix = get_matrix(&side);
    printf("Initial matrix:\n");
    __print_matrix(matrix, side);
    size_t op_len = side;
    Operation *ops = calloc(op_len, sizeof(Operation));
    size_t op_i = 0;
    if (straight_move(matrix, side, &ops, &op_len, &op_i) == 1) {
        printf("No non-zero element in column, exiting...\n");
        __print_matrix(matrix, side);
        return 1;
    }
    if (backwards_move(matrix, side, &ops, &op_len, &op_i) == 1) {
        printf("Zero element on diagonal in backward move, exiting...\n");
        __print_matrix(matrix, side);
        return 1;
    }
    printf("Result matrix:\n");
    __print_matrix(matrix, side);
    printf("Answers:\n");
    __print_answers(matrix, side);
    printf("Determinant: %lg\n", DETERMINANT);
    free(matrix[0]);
    free(matrix);
    op_len = op_i;
    ops = realloc(ops, op_len * sizeof(Operation));
    matrix = get_identity_matrix(side);
    int stat = apply_operations(matrix, side, ops, op_len);
    if (stat == SUBTRACT) {
        printf("No non-zero element in column, exiting...\n");
        __print_matrix(matrix, side);
        return 1;
    }
    if (stat == SUBTRACT_B) {
        printf("Zero element on diagonal in backward move, exiting...\n");
        __print_matrix(matrix, side);
        return 1;
    }
    printf("Inverse matrix:\n");
    __print_matrix(matrix, side);
    free(matrix[0]);
    free(matrix);
    free(ops);
    return 0;
}