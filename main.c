#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#define _MAIN_ELEM

double DETERMINANT = 1;

enum OPERATIONS
{
    NONE = 0,
    SWAP = 1,
    SUBTRACT = 2,
    SUBTRACT_B = 3
};

typedef struct operation
{
    int type;
    size_t which;
    size_t from;
} Operation;

void add_operation(Operation **ops, size_t *op_len, size_t *op_i, int type, size_t from, size_t which) {
    if (ops) {
        if (*op_i == *op_len) {
            *op_len *= 2;
            *ops = realloc(*ops, *op_len * sizeof(Operation));
        }
        (*ops)[*op_i].type = type;
        (*ops)[*op_i].which = which;
        (*ops)[*op_i].from = from;
        (*op_i)++;
    }
    return;
}

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
                matrix[i][j] = pow(q, i + j + 2) + 0.1 * ((int)(j - i));
            } else {
                matrix[i][j] = pow(q - 1, i + j + 2);
            }
        }
    }
    for (i = 0; i < side; i++) {
        matrix[i][side] = x * exp(x/(i+1)) * cos(x/(i+1));
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
        matrix[i] = matrix[0] + i * (*side+1);
    }
    if (!flag) {
        input_fill_matrix(matrix, *side);
    } else {
        formula_filled_matrix(matrix, *side);
    }
    return matrix;
}

void __print_matrix(double **matrix, size_t side) {
    size_t i, j;
    for (i = 0; i < side; i++) {
        printf("|");
        for (j = 0; j < side+1; ++j) {
            printf("%12lg|", matrix[i][j]);
        }
        printf("\n");
    }
    printf("#%-10s#\n", "");
    return;
}

int normalize_line(double **matrix, size_t side, size_t which) {
    size_t j, first_elem = (size_t)-1;
    for (j = 0; j < side; j++) {
        if (matrix[which][j]) {
            first_elem = j;
            break;
        }
    }
    if (first_elem == (size_t)-1) {
        return -1;
    }
    for (j = first_elem + 1; j < side+1; j++) {
        matrix[which][j] /= matrix[which][first_elem];
    }
    DETERMINANT *= matrix[which][first_elem];
    matrix[which][first_elem] = 1;
    return 0;
}

size_t find_fitting_row(double **matrix, size_t side, size_t column) {
    size_t j, res = (matrix[column][column] == 0 ? (size_t)-1 : column);
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

int subtract_line_from_line(double **matrix, size_t side, size_t from, size_t which) {
    size_t j, first_elem = (size_t)-1;
    for (j = 0; j < side; j++) {
        if (matrix[which][j]) {
            first_elem = j;
            break;
        }
    }
    if (first_elem == (size_t)-1) {
        return -1;
    }
    if (which != from) {
        double modifier = matrix[from][first_elem] / matrix[which][first_elem];
        for (j = first_elem + 1; j < side + 1; j++) {
            matrix[from][j] -= modifier * matrix[which][j];
        }
        matrix[from][first_elem] = 0;
    }
    return 0;
}

int subtract_line_from_line_b(double **matrix, size_t side, size_t from, size_t which) {
    if (!matrix[which][which]) {
        return 1;
    }
    matrix[from][side] -= matrix[which][side]*matrix[from][which];
    matrix[from][which] = 0;
    return 0;
}

int straight_move(double **matrix, size_t side, Operation **ops, size_t *op_len, size_t *op_i) {
    size_t column, row, s_row;
    for (column = 0; column < side; column++) {
        if ((s_row = find_fitting_row(matrix, side, column)) == (size_t)-1) {
            return 1;
        }
        for (row = column; row < side; row++) {
            if ((subtract_line_from_line(matrix, side, row, s_row) == -1)) {
                return 1;
            }
            if (ops) {
                add_operation(ops, op_len, op_i, SUBTRACT, row, s_row);
            }
        }
        swap_lines(matrix, side, s_row, column);
        if (ops) {
            add_operation(ops, op_len, op_i, SWAP, s_row, column);
        }
        normalize_line(matrix, side, column);
    }
    printf("Straight exit\n");
    return 0;
}

int backwards_move(double **matrix, size_t side, Operation **ops, size_t *op_len, size_t *op_i) {
    size_t s_row, row;
    for (s_row = side-1; s_row; s_row--) {
        for (row = s_row-1; row != (size_t)-1; row--) {
            if (subtract_line_from_line_b(matrix, side, row, s_row)) {
                return 1;
            }
            if (ops) {
                add_operation(ops, op_len, op_i, SUBTRACT_B, row, s_row);
            }
        }
    }
    return 0;
}

int main() {
#ifdef _MAIN_ELEM
    printf("Compiled with _MAIN_ELEM\n");
#else
    printf("Compiled without _MAIN_ELEM\n");
#endif
    size_t side;
    double **matrix = get_matrix(&side);
    __print_matrix(matrix, side);
    size_t op_len = side;
    Operation *ops = calloc(op_len, sizeof(Operation));
    size_t op_i = 0;
    if (straight_move(matrix, side, &ops, &op_len, &op_i) == 1) {
        printf("No non-zero element in column, exiting...\n");
        __print_matrix(matrix, side);
        return 1;
    }
    __print_matrix(matrix, side);
    if (backwards_move(matrix, side, &ops, &op_len, &op_i) == 1) {
        printf("Zero element on diagonal in backward move, exiting...\n");
        __print_matrix(matrix, side);
        return 1;
    }
    __print_matrix(matrix, side);
    printf("Determinant: %lg\n", DETERMINANT);
    free(ops);
    free(matrix[0]);
    free(matrix);
    return 0;
}