#include <cstdio>
#include <stdlib.h>
#include <string.h>
#include <cmath>

static inline void copy(const double *source, double *dest, const int v, const int h){
    memcpy(dest, source, v * h * sizeof(double));
}

#define max(x,y) (x > y ? x : y)
#define min(a, b) (a < b ? a : b)
#define A(i,j) (A + i*n*m + j*av*m)
#define a(p,q) (A(i,j) + p*ah + q)
#define eps (1e-15) // нормы для невязки правой части и решения

void f0(double *const A, const int n, const int m, const char* const filename, int *errno)
{
    int i, j, p, q, av, ah;
    const int k = n / m;
    const int l = n - k * m;
    int counter = 0;
    FILE* file_input = fopen(filename, "r");
    if (file_input == NULL) {
        *errno = 3;
        return;
    }

    for (i = 0; i * m < n; i++) {
        av = i < k ? m : l;
        for (p = 0; p < av; p++) {
            for (j = 0; j * m < n; j++) {
            ah = j < k ? m : l;
                for (q = 0; q < ah; q++) {
                    if (fscanf(file_input, "%lf", a(p, q)) == 1)
                        counter++;
                    else{
                        break;
                    }
                }
            }
        }
    }
    if (counter != n * n) *errno = 4;
    fclose(file_input);
}
void f1(double *const A, const int n, const int m) {
    int i, j, p, q, av, ah;
    const int k = n / m;
    const int l = n - k * m;

    for (i = 0; i * m < n; i++) {
        av = i < k ? m : l;
        for (j = 0; j * m < n; j++) {
            ah = j < k ? m : l;
            for (p = 0; p < av; p++) {
                for (q = 0; q < ah; q++) {
                    *a(p, q) = n - max(i * m + p + 1, j * m + q + 1) + 1; // i*m + p + 1 это номер строки этого элемента во всей исходной матрице (+1 потому что нумерация с 1)
                }
            }
        }
    }
}
void f2(double *const A, const int n, const int m) {
    int i, j, p, q, av, ah;
    const int k = n / m;
    const int l = n - k * m;

    for (i = 0; i * m < n; i++) {
        av = i < k ? m : l;
        for (j = 0; j * m < n; j++) {
            ah = j < k ? m : l;
            for (p = 0; p < av; p++) {
                for (q = 0; q < ah; q++) {
                    *a(p, q) = max(i * m + p + 1, j * m + q + 1);
                }
            }
        }
    }
}
void f3(double *const A, const int n, const int m) {
    int i, j, p, q, av, ah;
    const int k = n / m;
    const int l = n - k * m;

    for (i = 0; i * m < n; i++) {
        av = i < k ? m : l;
        for (j = 0; j * m < n; j++) {
            ah = j < k ? m : l;
            for (p = 0; p < av; p++) {
                for (q = 0; q < ah; q++) {
                    *a(p, q) = abs(i * m + p - (j * m + q));
                }
            }
        }
    }
}
void f4(double* const A, const int n, const int m) {
    int i, j, p, q, av, ah;
    const int k = n / m;
    const int l = n - k * m;

    for (i = 0; i * m < n; i++)
    {
        av = i < k ? m : l;
        for (j = 0; j * m < n; j++)
        {
            ah = j < k ? m : l;
            for (p = 0; p < av; p++)
            {
                for (q = 0; q < ah; q++) *a(p, q) = fabs(1. / (i * m + p + j * m + q + 1));
            }
        }
    }
}

void fill(double *const A, const int n, const int m, const int s, const char *const filename, int *const errno)
{
    switch (s) {
    case 0: f0(A, n, m, filename, errno); break;
    case 1: f1(A, n, m); break;
    case 2: f2(A, n, m); break;
    case 3: f3(A, n, m); break;
    case 4: f4(A, n, m); break;
    default: break;
    }
}

void fill_right_part(const double *const A, double *const B, const int n, const int m) // вероятно это не нужно
{
    int i, j, p, q, av, ah;
    const int k = n / m;
    const int l = n - k * m;
    int c = 0;

    for (i = 0; i < n; i++) B[i] = 0; // вероятно это необязательно
    
    for (i = 0; i * m < n; i++) { // цикл по строкам
        av = i < k ? m : l;
        for (p = 0; p < av; p++) {
            B[i * m + p] = 0;
            c = 0;
            for (j = 0; j * m < n; j++){
                // здесь было c = 0;
                ah = j < k ? m : l;
                for (q = 0; q < ah; q++, c++){
                    //printf("%d  ", c);
                    /* if (c % 2 == 0) */ B[i * m + p] += *a(p, q) * ((c + 1) % 2); // напоминание: четность другая, потому что нумерация в матрице и в массиве идет с 1 и 0 соответственно
                }
            }
        }
    }
}

double full_norm(const double *const A, const int n, const int m){
    const int k = n / m;
    const int l = n - k * m;
    int i, j, p, q, av, ah;
    double max = 0., current = 0.;

    for (i = 0; i * m < n; i++) {
        av = i < k ? m : l;
        for (p = 0; p < av; p++) {
            current = 0.;
            for (j = 0; j * m < n; j++){
                ah = j < k ? m : l;
                for (q = 0; q < ah; q++){
                    current += fabs(*a(p, q));
                }
            }
            if (current > max) max = current;
        }
    }
    return max;
}

double st_norm(const double *const A, const int n){
    int i, j;
    double max = 0., current = 0.;
    for (i = 0; i < n; i++) {
        current = 0.;
        for (j = 0; j < n; j++) {
            current += fabs(A[j * n + i]);
        }
        if (current > max)
            max = current;
    }
    return max;
}
// error
int error(const int error) {
    switch (error) {
    case 1: {
        printf("Некорректный ввод команды\n");
        break;
    }
    case 2: {
        printf("Ошибка выделения памяти\n");
        break;
    }
    case 3: {
        printf("Ошибка чтения файла\n");
        break;
    }
    case 4: {
        printf("Некорректные данные или недостаточно данных\n");
        break;
    }
    case 5: {
        printf("Алгоритм неприменим\n");
        break;
    }
    default:
        return 0;
    }
    return error;
}

double res1(double *A, double *B, double *X, const int n, const int m) {
    int i, j, q, r;
    int av, ah;
    double *pa, *pi, *pj;
    const int k = n / m;
    const int l = n - k * m;
    double norm = 0.;
    double sum = 0.;

    for (i = 0; i < n; i++) {
        norm += fabs(B[i]);
    }
    for (j = 0; j * m < n; j++) {
        ah = (j < k) ? m : l;
        pj = X + j * m;
        for (i = 0; i * m < n; i++) {
            pi = B + i * m;
            av = (i < k) ? m : l;
            pa = A + i * m * n + j * av * m;
            for (q = 0; q < av; q++) {
                sum = 0.;
                for (r = 0; r < ah; r++) {
                    sum += pa[q * ah + r] * pj[r];
                }
                pi[q] -= sum;
            }
        }
    }
    sum = 0.;
    if (norm > eps) {
        for (i = 0; i < n; i++) {
            sum += fabs(B[i]);
        }
        sum /= norm;
    }

    return sum;
}

double res2(double *B, double *X, const int n) { // B это вектор в знаменателе
    int i;
    double norm1 = 0., norm2 = 0.;

    for (i = 0; i < n; i++) {
        //norm1 += ((i + 1) % 2);
        norm1 += fabs(B[i]);
    }
    //if (norm1 > eps) {
    for (i = 0; i < n; i++) {
        //norm2 += fabs(B[i] - ((i + 1) % 2));
        norm2 += fabs(B[i] - X[i]);
    }
    norm2 /= norm1;
    //}
    return norm2;
}

void print_matrix(const double *A, const int h, const int v, const int m, const int r) { // вывод матрицы
    int nv = min(v, r);
    int nh = min(h, r); // если матрица будет меньше, чем размер вывода, то это условие будет существенно

    int mv, mh;

    int prn_val_h = r;
    int prn_val_v = r;

    int i, j, p, q, av, ah;
    int kh = h / m, kv = v / m;
    int lh = nh - kh * m, lv = nv - kv * m;

    for (i = 0; i * m < nh; i++) {
        av = i < kh ? m : lh;
        mv = min(av, prn_val_v);
        for (p = 0; p < mv; p++) {
            prn_val_h = r;
            for (j = 0; j * m < nv; j++) {
                ah = j < kv ? m : lv;
                mh = min(ah, prn_val_h);
                for (q = 0; q < mh; q++) printf(" %10.3e", *(A + i*h*m + j*av*m  + p*ah + q));
                prn_val_h -= ah;
            }
            printf("\n");
        }
        prn_val_v -= av;
    }
    printf("\n");
}

// c = a * b
int multiply(const double* const a, const int av, const int ah, const double* const b, const int bv, const int bh, double* const c) {
    int r = 0, t = 0, q, temp_r, temp_t;
    // для вышеописанной логики
    double c00, c01, c02, c10, c11, c12, c20, c21, c22;
    if (ah != bv)
        return -1;
    // Зануляем c
    for (r = 0; r < av; r++)
        for (t = 0; t < bh; t++)
            c[r * bh + t] = 0.;

    for (r = 0; r + 2 < av; r += 3)
        for (t = 0; t + 2 < bh; t += 3) {
            c00 = 0., c01 = 0., c02 = 0.;
            c10 = 0., c11 = 0., c12 = 0.;
            c20 = 0., c21 = 0., c22 = 0.;
            for (q = 0; q < ah; q++) {
                c00 += a[(r + 0) * ah + q] * b[q * bh + (t + 0)];
                c01 += a[(r + 0) * ah + q] * b[q * bh + (t + 1)];
                c02 += a[(r + 0) * ah + q] * b[q * bh + (t + 2)];
                c10 += a[(r + 1) * ah + q] * b[q * bh + (t + 0)];
                c11 += a[(r + 1) * ah + q] * b[q * bh + (t + 1)];
                c12 += a[(r + 1) * ah + q] * b[q * bh + (t + 2)];
                c20 += a[(r + 2) * ah + q] * b[q * bh + (t + 0)];
                c21 += a[(r + 2) * ah + q] * b[q * bh + (t + 1)];
                c22 += a[(r + 2) * ah + q] * b[q * bh + (t + 2)];
            }
            c[(r + 0) * bh + (t + 0)] += c00;
            c[(r + 0) * bh + (t + 1)] += c01;
            c[(r + 0) * bh + (t + 2)] += c02;
            c[(r + 1) * bh + (t + 0)] += c10;
            c[(r + 1) * bh + (t + 1)] += c11;
            c[(r + 1) * bh + (t + 2)] += c12;
            c[(r + 2) * bh + (t + 0)] += c20;
            c[(r + 2) * bh + (t + 1)] += c21;
            c[(r + 2) * bh + (t + 2)] += c22;
        }
    temp_t = t;
    temp_r = r;
    // если не получилось разделить четно, то
    // повторяем процесс для последнего столбца и строчки
    // так, как делали раньше
    if ((av - temp_r) == 2) {
        for (t = 0; t + 1 < bh; t += 2) {
            c00 = 0., c01 = 0.;
            c10 = 0., c11 = 0.;
            for (q = 0; q < ah; q++) {
                c00 += a[(temp_r + 0) * ah + q] * b[q * bh + (t + 0)];
                c01 += a[(temp_r + 0) * ah + q] * b[q * bh + (t + 1)];
                c10 += a[(temp_r + 1) * ah + q] * b[q * bh + (t + 0)];
                c11 += a[(temp_r + 1) * ah + q] * b[q * bh + (t + 1)];
            }
            c[(temp_r + 0) * bh + (t + 0)] = c00;
            c[(temp_r + 0) * bh + (t + 1)] = c01;
            c[(temp_r + 1) * bh + (t + 0)] = c10;
            c[(temp_r + 1) * bh + (t + 1)] = c11;
        }
        if (t < bh) {
            for (r = temp_r; r < av; r++)
                for (t = bh & ~1; t < bh; t++) {
                    c00 = 0.;
                    for (q = 0; q < ah; q++) {
                        c00 += a[(r + 0) * ah + q] * b[q * bh + (t + 0)];
                    }
                    c[(r + 0) * bh + (t + 0)] = c00;
                }
        }
    }

    if ((bh - temp_t) == 2) {
        for (r = 0; r + 1 < av; r += 2) {
            c00 = 0., c01 = 0.;
            c10 = 0., c11 = 0.;
            for (q = 0; q < ah; q++) {
                c00 += a[(r + 0) * ah + q] * b[q * bh + (temp_t + 0)];
                c01 += a[(r + 0) * ah + q] * b[q * bh + (temp_t + 1)];
                c10 += a[(r + 1) * ah + q] * b[q * bh + (temp_t + 0)];
                c11 += a[(r + 1) * ah + q] * b[q * bh + (temp_t + 1)];
            }
            c[(r + 0) * bh + (temp_t + 0)] = c00;
            c[(r + 0) * bh + (temp_t + 1)] = c01;
            c[(r + 1) * bh + (temp_t + 0)] = c10;
            c[(r + 1) * bh + (temp_t + 1)] = c11;
        }
        for (r = av & ~1; r < av; r++) {
            for (t = temp_t; t < bh; t++) {
                c00 = 0.;
                for (q = 0; q < ah; q++) {
                    c00 += a[(r + 0) * ah + q] * b[q * bh + (t + 0)];
                }
                c[(r + 0) * bh + (t + 0)] = c00;
            }
        }
    }

    if (av - temp_r == 1) {
        for (t = 0; t < bh; t++) {
            c00 = 0.;
            for (q = 0; q < ah; q++)
                c00 += a[temp_r * ah + q] * b[q * bh + t];
            c[temp_r * bh + t] = c00;
        }
    }
    if (bh - temp_t == 1) {
        for (r = 0; r < av; r++) {
            c00 = 0.;
            for (q = 0; q < ah; q++)
                c00 += a[r * ah + q] * b[q * bh + temp_t];
            c[r * bh + temp_t] = c00;
        }
    }
    return 0;
}

// c -= a * b
// эта функция неверна, что-то в ней необходимо переписать

// c -= b
int extract(double *const c, const int cv, const int ch, const double *const b, const int bv, const int bh){ // вторая версия параллельная, но это не нужно, так как это для блоков, которые целиком умещаются в оперативную память
    if (cv != bv || ch != bh) return -1;
    int i, j;
    for(i = 0; i < cv; i++)
        for(j = 0; j < ch; j++) c[i * ch + j] -= b[i * ch + j];
    return 0;
}

#undef A

void swap(double* const lhs, double* const rhs) {
    double temp = *lhs;
    *lhs = *rhs;
    *rhs = temp;
}

#define A(i, j) A[i * n + j]
#define E(i, j) A_inversed[i * n + j]

int jord_inverse(double *A, double *A_inversed, const int n, double ERROR) { // ERROR тут фиктивен, надо его убрать
    int i, j; // итераторы
    double max = 0.; // максимальный элемент по столбцу
    int m = 0; // строчка и столбец с максимальным элементом
    double c = 0.; // коэффициент
    
    for (j = 0; j < n; j++) { // идем по столбам
        max = fabs(A(j, j));
        m = j;

        for (i = j; i < n; i++) {
            if (fabs(A(i, j)) > max) {
                m = i;
                max   = fabs(A(i, j));
            }
        } // нашли максимальный элемент и строку в которой он находится

        if (max <= ERROR) { // ошибка нужна чтобы нормально найти максимальный элемент
            return -1;
        }

        if (m != j) {
            for (i = 0; i < j; i++) { // меняем строки местами
                swap(&(E(j, i)), &(E(m, i))); // элемент на дигонали мы не меняем так как там все равно будет 1 в (j,j), а все остальные 0
            }
            for (i = j; i < n; i++) { // менять их в матрице А не надо так как там нули
                swap(&(A(j, i)), &(A(m, i)));
                swap(&(E(j, i)), &(E(m, i)));
            }
        }
        //c = 1. / A(j, j);
        c = A(j,j); // этот элемент на диагонали не равен 0
        //printf("%lf    ", c);
        //if (fabs(c - 1) > ERROR) {
            for (m = 0; m < j + 1; m++) E(j, m) /= c; // m больше не нужен так как строки уже переставили
            A(j, j) = 1.;
            for (m = j + 1; m < n; m++) {
                A(j, m) /= c;
                E(j, m) /= c;
            }
        //}

        // цикл ниже вычитает обнуляет почти весь столб
        for (i = 0; i < n; i++) { // все элементы в столбе j кроме одного равны 0
            if (i != j){
                c = A(i, j);
                //if (fabs(c) > ERROR) {
                    for (m = 0; m < j + 1; m++) {
                        E(i, m) -= c * E(j, m);
                    }
                    // A(i, j) = 0.;
                    for (m = j + 1; m < n; m++) {
                        A(i, m) -= c * A(j, m);
                        E(i, m) -= c * E(j, m);
                    }
                //}
            }
        }
    } // теперь в Rev находится обратная матрица решенная Жорданом
    return 0;
}

int solve(int n, int m, double *A, double *B, double *X, double *V1, double *V2, double *V3) { /////////////////////////////////////////////////////////////////////////////////////////
    int i = 0, j = 0, r = 0, q = 0;                     // итераторы
    double *pa = NULL, *pi = NULL, *pj = NULL;          // вспомогательные указатели
    int av = 0, ah = 0;                                 // размер текущего блока av * ah
    const int k = n / m;                                // количество блоков
    const int l = n - k * m;                            // крайний блок имеет сторону ...
    double ERROR = (full_norm(A, n, m) * eps);          // примерная допустимая погрешность
    double min = 0.;                                    // минимальная норма обратной матрицы (по столбцам)
    int  min_i = 0;                                     // номер строки с минимальной нормой обратной матрицы в столбце
    int c = 0;                                          // счетчик необратимых матриц в столбце
    double current = 0.;                                // текущая норма

    for (j = 0; j * m < n; j++) { // в этом цикле смотрим квадратные блоки порядка m. 
        av = ah = j < k ? m : l;
        pa = A + j * n * m + j * av * m;

        // A_{j, j} --> V1
        // V_min = (A_{j, j})^(-1)
        // min = ||V_min||
        for (r = 0; r < av; r++) { // тут работаем только с первым блоком (делаем первый шаг, чтобы записать указатели V1, V2, V3)
            for (q = 0; q < ah; q++) {
                V1[r * ah + q] = pa[r * ah + q];
                V3[r * ah + q] = (r == q);
            }
        }
        if (jord_inverse(V1, V3, av, ERROR) == 0) { // ищем матрицу с наименьшей нормой обратного
            min   = st_norm(V3, av);
            min_i = j;
        } else {
            min = 0.;
            c++;
        }

        // i = j + 1, ..., k - 1
        // A_{i, j} --> V1, V2 = (A_{i, j})^(-1)
        // if ||V2|| < min V3 = V2, min = ||V2||
        for (i = j + 1; i * m + l < n; i++) {
            pi = A + i * n * m + j * av * m; // это блок (i, j)

            for (r = 0; r < av; r++) { // вносим pi в указатель V1
                for (q = 0; q < ah; q++) {
                    V1[r * ah + q] = pi[r * ah + q];
                    V2[r * ah + q] = (r == q);
                }
            }
            if (jord_inverse(V1, V2, av, ERROR) == 0) { // если существует обратная у V1
                current = st_norm(V2, av);
                if (current < min) {
                    pj = V2;
                    V2 = V3;
                    V3 = pj;
                    min   = current;
                    min_i = i;
                }
            } else {
                c++;
            }
        }
        
        if (c == k - j + int(av == l)) { // если ни один блок не обратим, то алгоритм неприменим
            return -1;
        }

        if (min_i != j) { // меняем местами строки блоков
            for (i = 0; i * m < n; i++) {
                q  = (i < k) ? m : l; // мы помним что нумерация блоков идет с 0, поэтому столб k содержит неквадратные блоки
                pi = A + min_i * n * m + i * q * m;
                pj = A + j * n * m + i * q * m; // here is my fucking m instead q

                for (r = 0; r < m; r++) { // меняем местами блоки в i столбе?
                    for (c = 0; c < q; c++) {
                        current       = pi[r * q + c];
                        pi[r * q + c] = pj[r * q + c];
                        pj[r * q + c] = current;
                    }
                }
            }
            pi = B + min_i * m;
            pj = B + j * m;

            for (r = 0; r < m; r++) { // меняем местами блоки свободного вектора
                current = pi[r];
                pi[r]   = pj[r];
                pj[r]   = current;
            }
        }

        // A_{j, j} = E, V_3 * (A_{j, j+1},...,A_{j,k+1},B_{j})
        // E_matrix(pa, av);
        for (i = j + 1; i * m < n; i++) { // for (i = 0; i * m < n; i++) { 
            r = (i < k) ? m : l;
            pi = A + j * n * m + i * av * m;
            copy(pi, V1, av, r);
            multiply(V3, av, av, V1, av, r, V2);
            copy(V2, pi, av, r);
        }
        pi = B + j * m;
        copy(pi, V1, av, 1);
        multiply(V3, av, av, V1, av, 1, V2);
        copy(V2, pi, av, 1);
/*         print_matrix(A,n,n,m,n);
        print_matrix(B,1,n,m,n); */

        // A_{ i, c } = A_{ i, c } - A_{ i, j } x A_{ j, c }
        //      pa    =     pa     -     pi     x     pj
        //for (i = 0; i < k; i++){
        for (i = 0; i * m < n; i++) { // идем по строчкам вниз
            if (i != j){
                /* if (i == k) printf("%d\n", i); */
                q = (i < k) ? m : l;
                pi = A + i * n * m + j * q * m;
                copy(pi, V1, q, av);

                for (c = j + 1; c * m < n; c++) { // каждую умножаем и вычитаем с подходящим коэффицентом
                    ah = (c < k) ? m : l;
                    pa = A + i * n * m + c * q * m;
                    pj = A + j * n * m + c * av * m;

                    copy(pj, V2, av, ah);
                    multiply(V1, q, av, V2, av, ah, V3); // было изменено здесь, вместо av, m
                    extract(pa, q,ah,V3,q,ah);
                    // multiply_and_extract(V1, q, av, V2, av, ah, pa);
                }
                pa = B + i * m;
                pj = B + j * m;
                copy(pj, V2, av, 1);
                multiply(V1, q, av, V2, av, 1, V3); // было изменено, здесь вместо av, m
                extract(pa, q,1,V3,q,1);
                // multiply_and_extract(V1, q, av, V2, av, 1, pa);
            }
        }
        c = 0;
    }
    //print_matrix(A,n,n,m,n);
    copy(B, X, n, 1); // в этой строке можно убрать, ибо решение записано в B
    //print_matrix(B,1,n,m,n);
    return 0;
}

#undef eps