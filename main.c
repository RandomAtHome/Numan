#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double DETERMINANT = 1;

double **get_matrix(size_t side) {
    double **matrix = calloc(side, sizeof(double *));
    matrix[0] = calloc(side * (side + 1), sizeof(double));
    size_t i;
    for (i = 1; i < side; i++) {
        matrix[i] = matrix[0] + i * (side+1);
    }
    return matrix;
}

void fill_matrix(double **matrix, size_t side) {
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

void __print_matrix(double **matrix, size_t side) {
    size_t i, j;
    for (i = 0; i < side; i++) {
        for (j = 0; j < side+1; ++j) {
            printf("%lf\t", matrix[i][j]);
        }
        printf("\n");
    }
    printf("#%-10s#\n", "");
    return;
}

double fabs(double number) {
    return number * (number > 0 ? 1 : -1);
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
    normalize_line(matrix, side, which);
    return 0;
}

int straight_move(double **matrix, size_t side) {
    size_t column, row, s_row;
    for (column = 0; column < side; column++) {
        if ((s_row = find_fitting_row(matrix, side, column)) == (size_t)-1) {
            return 1;
        }
        for (row = column; row < side; row++) {
            if ((subtract_line_from_line(matrix, side, row, s_row) == -1)) {
                return 1;
            }
        }
        if (column != s_row) {
            swap_lines(matrix, side, s_row, column);
        }
    }
    return 0;
}

int backwards_move(double **matrix, size_t side) {
    size_t s_row, row;
    for (s_row = side-1; s_row; s_row--) {
        if (!matrix[s_row][s_row]) {
            return 1;
        }
        for (row = s_row-1; row != (size_t)-1; row--) {
            matrix[row][side] -= matrix[s_row][side]*matrix[row][s_row];
            matrix[row][s_row] = 0;
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
    printf("Enter n - side of matrix:\n");
    scanf("%zu", &side);
    double **matrix = get_matrix(side);
    fill_matrix(matrix, side);
    __print_matrix(matrix, side);
    if (straight_move(matrix, side) == 1) {
        printf("No non-zero element in column, exiting...\n");
        return 1;
    }
    __print_matrix(matrix, side);
    if (backwards_move(matrix, side) == 1) {
        printf("Zero element on diagonal in backward move, exiting...\n");
        return 1;
    }
    __print_matrix(matrix, side);
    printf("%lf\n", DETERMINANT);
    return 0;
}