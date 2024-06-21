#include "matrix.h"
Matrix createMatrix(int rows, int cols) {
    Matrix matrix;
    matrix.rows = rows;
    matrix.cols = cols;
    matrix.data.resize(rows, vector<double>(cols));
    return matrix;
}
void fillMatrix(Matrix& matrix) {
    cout << "Введите элементы матрицы:\n";
    for (int i = 0; i < matrix.rows; ++i) {
        for (int j = 0; j < matrix.cols; ++j) {
            cout << "Введите элемент [" << i + 1 << "][" << j + 1 << "]: ";
            double element = getValidInput("");
            matrix.data[i][j] = element;
        }
    }
}

double getValidInput(const string& prompt) {
    string input;
    double value;
    cout << prompt;
    while (true) {
        getline(cin, input);

        size_t start = input.find_first_not_of(" ");
        if (start == string::npos) {
            cout << "Ошибка: Неверный формат числа. Введите снова.\n";
            continue;
        }

        bool valid = true;
        bool pointFound = false; 
        bool signFound = false; 

        for (size_t i = start; i < input.length(); ++i) {
            char c = input[i];

            if (!isdigit(c)) {
                if (c == '.' && !pointFound) { 
                    pointFound = true;
                }
                else if ((c == '+' || c == '-') && i == start && !signFound) {
                    signFound = true;
                }
                else {
                    valid = false;
                    break;
                }
            }
        }

        if (valid) {
            try {
                value = stod(input);
                break;
            }
            catch (const invalid_argument& e) {
                cout << "Ошибка: Неверный формат числа. Введите снова.\n";
            }
            catch (const out_of_range& e) {
                cout << "Ошибка: Число слишком большое. Введите снова.\n";
            }
        }
        else {
            cout << "Ошибка: Неверный формат числа. Введите снова.\n";
        }
    }
    return value;
}

double calculateDeterminant(const Matrix& matrix) {
    if (matrix.rows != matrix.cols) {
        throw invalid_argument("Матрица должна быть квадратной");
    }

    int n = matrix.rows;
    if (n == 1) {
        return matrix.data[0][0];
    }
    else if (n == 2) {
        return matrix.data[0][0] * matrix.data[1][1] - matrix.data[0][1] * matrix.data[1][0];
    }
    else {
        double det = 0;
        for (int i = 0; i < n; ++i) {
            Matrix submatrix = createMatrix(n - 1, n - 1);
            int subi = 0;
            for (int j = 1; j < n; ++j) {
                int subj = 0;
                for (int k = 0; k < n; ++k) {
                    if (k != i) {
                        submatrix.data[subi][subj] = matrix.data[j][k];
                        ++subj;
                    }
                }
                ++subi;
            }
            det += pow(-1, i) * matrix.data[0][i] * calculateDeterminant(submatrix);
        }
        return det;
    }
}

Matrix getInverseMatrix(const Matrix& matrix) {
    int n = matrix.rows;
    Matrix augmented = createMatrix(n, 2 * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented.data[i][j] = matrix.data[i][j];
        }
        augmented.data[i][i + n] = 1;
    }
    for (int i = 0; i < n; ++i) {
        double glav = augmented.data[i][i];
        if (fabs(glav) < 1e-9) {
            cout << "Матрица вырождена, обратная матрица не существует." << endl;
            return Matrix(); 
        }
        for (int j = 0; j < 2 * n; ++j) {
            augmented.data[i][j] /=glav;
        }
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = augmented.data[k][i];
                for (int j = 0; j < 2 * n; ++j) {
                    augmented.data[k][j] -= factor * augmented.data[i][j];
                }
            }
        }
    }

    Matrix inverse = createMatrix(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse.data[i][j] = augmented.data[i][j + n];
        }
    }

    return inverse;
}

vector<double> solveCramer(const Matrix& A, const vector<double>& b) {
    int n = A.rows;
    double detA = calculateDeterminant(A);

    if (detA == 0) {
        cerr << "Ошибка: Определитель матрицы равен нулю. Решение невозможно." << endl;
        return {};
    }

    vector<double> x(n);
    for (int i = 0; i < n; ++i) {
        Matrix tempA = A;
        for (int j = 0; j < n; ++j) {
            tempA.data[j][i] = b[j];
        }
        x[i] = calculateDeterminant(tempA) / detA;
    }
    return x;
}

vector<double> solveInverseMatrix(const Matrix& A, const vector<double>& b) {
    int n = A.rows;
    Matrix A_inv = getInverseMatrix(A);

    if (A_inv.data.empty()) {
        cerr << "Ошибка: Матрица необратима. Решение невозможно." << endl;
        return {};
    }

    vector<double> x(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            x[i] += A_inv.data[i][j] * b[j];
        }
    }
    return x;
}
void showMenu() {
    int choice;
    Matrix A;
    vector<double> b, x;
    ofstream file1("baza.txt", ios::app);

    do {
        cout << "\nМеню:\n";
        cout << "1. Решить СЛАУ методом Крамера\n";
        cout << "2. Решить СЛАУ методом обратной матрицы\n";
        cout << "3. Работа с файлом\n";
        cout << "0. Выход\n";
        cout << "Выберите действие: ";
        cin >> choice;
        cin.ignore(numeric_limits<streamsize>::max(), '\n');

        switch (choice) {
        case 1: {
            int n;
            cout << "Введите размер матрицы: ";
            string input;
            while (true) {
                getline(cin, input);
                bool valid = true;
                for (char c : input) {
                    if (!isdigit(c)) {
                        valid = false;
                        break;
                    }
                }
                if (valid) {
                    n = stoi(input); 
                    if (n > 0) {
                        break;
                    }
                }
                cout << "Ошибка: Неверный размер матрицы. Введите положительное целое число: ";
            }
            A = createMatrix(n, n);
            fillMatrix(A);
            b.resize(n);
            cout << "Введите свободные члены:\n";
            for (int i = 0; i < n; ++i) {
                b[i] = getValidInput("Введите свободный член " + to_string(i + 1) + ": ");
            }
            auto start = chrono::high_resolution_clock::now();

            x = solveCramer(A, b);
            if (!x.empty()) {
                cout << "Решение:\n";
                for (int i = 0; i < n; ++i) {
                    cout << "x[" << i+1 << "] = " << x[i] << endl;
                }
            }
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
            cout << "Время выполнения: " << duration << " микросекунд" << endl;
            file1 << "Решение методом Крамера:\n";
            for (int i = 0; i < n; ++i) {
                file1 << "x[" << i + 1 << "] = " << x[i] << endl;
            }
            file1 << "Время выполнения: " << duration << " микросекунд\n\n";
            break;  
        }
        case 2: {
            int n;
            cout << "Введите размер матрицы: ";
            string input;
            while (true) {
                getline(cin, input);
                bool valid = true;
                for (char c : input) {
                    if (!isdigit(c)) {
                        valid = false;
                        break;
                    }
                }
                if (valid) {
                    n = stoi(input);
                    if (n > 0) {
                        break;
                    }
                }
                cout << "Ошибка: Неверный размер матрицы. Введите положительное целое число: ";
            }
            A = createMatrix(n, n);
            fillMatrix(A);
            b.resize(n);
            for (int i = 0; i < n; ++i) {
                b[i] = getValidInput("Введите свободный член " + to_string(i + 1) + ": ");
            }
            auto start = chrono::high_resolution_clock::now();
            x = solveInverseMatrix(A, b);
            if (!x.empty()) {
                cout << "Решение:\n";
                for (int i = 0; i < n; ++i) {
                    cout << "x[" << i+1 << "] = " << x[i] << endl;
                }
            }
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
            cout << "Время выполнения: " << duration << " микросекунд" << endl;
            file1 << "Решение методом обратной матрицы:\n";
            for (int i = 0; i < n; ++i) {
                file1 << "x[" << i + 1 << "] = " << x[i] << endl;
            }
            file1 << "Время выполнения: " << duration << " микросекунд\n\n";

            break;
        }
        case 3: {
            string filename;
            cout << "Введите имя файла с примером: ";
            getline(cin, filename);

            ifstream file(filename);
            if (file.is_open()) {
                int n;
                file >> n; 
                A = createMatrix(n, n);
                b.resize(n);
                for (int i = 0; i < n; ++i) {
                    for (int j = 0; j < n; ++j) {
                        file >> A.data[i][j];
                    }
                }

                for (int i = 0; i < n; ++i) {
                    file >> b[i];
                }
                file.close();
                cout << "Выберите метод решения:\n";
                cout << "1. Метод Крамера\n";
                cout << "2. Метод обратной матрицы\n";
                int methodChoice;
                cin >> methodChoice;
                cin.ignore(numeric_limits<streamsize>::max(), '\n');

                if (methodChoice == 1) {
                    auto start = chrono::high_resolution_clock::now();
                    x = solveCramer(A, b);
                    if (!x.empty()) {
                        cout << "Решение:\n";
                        for (int i = 0; i < n; ++i) {
                            cout << "x[" << i + 1 << "] = " << x[i] << endl;
                        }

                    }
                    auto end = chrono::high_resolution_clock::now();
                    auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
                    cout << "Время выполнения: " << duration << " микросекунд" << endl;
                    file1 << "Решение методом Крамера:\n";
                    for (int i = 0; i < n; ++i) {
                        file1 << "x[" << i + 1 << "] = " << x[i] << endl;
                    }
                    file1 << "Время выполнения: " << duration << " микросекунд\n\n";
                }
                else if (methodChoice == 2) {
                    auto start = chrono::high_resolution_clock::now();
                    x = solveInverseMatrix(A, b);
                    if (!x.empty()) {
                        cout << "Решение:\n";
                        for (int i = 0; i < n; ++i) {
                            cout << "x[" << i + 1 << "] = " << x[i] << endl;
                        }
                        auto end = chrono::high_resolution_clock::now();
                        auto duration = chrono::duration_cast<chrono::microseconds>(end - start).count();
                        cout << "Время выполнения: " << duration << " микросекунд" << endl;
                        file1 << "Решение методом обратной матрицы:\n";
                        for (int i = 0; i < n; ++i) {
                            file1 << "x[" << i + 1 << "] = " << x[i] << endl;
                        }
                        file1 << "Время выполнения: " << duration << " микросекунд\n\n";
                    }
                    break;

                }
                else {
                    cout << "Неверный выбор метода.\n";
                }
            }
            else {
                cout << "Ошибка: Не удалось открыть файл " << filename << endl;
            }
            break;
        }
        case 0:
            cout << "Программа завершена.\n";
            break;
        default:
            cout << "Неверный выбор. Попробуйте снова.\n";
        }
    } while (choice != 0);
}