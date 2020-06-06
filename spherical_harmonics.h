#define PI 3.14159265359

unsigned long int fact(int n) {

  unsigned long int result = 1;
  while (n) {
    result = n * result;
    --n;
  }
  return result;
}

double K_lm(int l, int m) {
  return sqrt((2 * l + 1) * fact(l - abs(m)) / (4 * PI * fact(l + abs(m))));
}

double P_lm(int l, int m, double theta) {
  double x = cos(theta);

  if (l == 0 && m == 0)
    return 1;
  else if (l == 1 && m == 0)
    return x;
  else if (l == 1 && m == 1 || m == -1)
    return -1.0 * sqrt(1 - x * x);
  else if (l == 2 && m == 0)
    return 0.5 * (3.0 * x * x - 1);
  else if (l == 2 && m == 1 || m == -1)
    return -3.0 * x * sqrt(1 - x * x);
  else if (l == 2 && m == 2 || m == -2)
    return 3.0 * (1 - x * x);
  else if (l == 3 && m == 0)
    return 0.5 * x * (5.0 * x * x - 3);
  else if (l == 3 && m == 1 || m == -1)
    return 1.5 * (1 - 5.0 * x * x) * sqrt(1 - x * x);
  else if (l == 3 && m == 2 || m == -2)
    return 15.0 * x * (1 - x * x);
  else if (l == 3 && m == 3 || m == -3)
    return -15.0 * pow((1 - x * x), 1.5);
  else if (l == 4 && m == 0)
    return (1.0 / 8.0) * (35.0 * x * x * x * x - 30.0 * x * x + 3.0);
  else if (l == 4 && m == 1 || m == -1)
    return (5 / 2) * x * (3 - 7 * x * x) * sqrt(1 - x * x);
  else if (l == 4 && m == 2 || m == -2)
    return (15 / 2) * (7.0 * x * x - 1) * (1 - x * x);
  else if (l == 4 && m == 3 || m == -3)
    return -105.0 * x * pow((1 - x * x), 1.5);
  else if (l == 4 && m == 4 || m == -4)
    return 105 * pow((1 - x * x), 2.0);

}

double y_lm(double phi, double theta, int l, int m) {

  if (m > 0)
    return sqrt(2.0) * K_lm(l, m) * cos(m * phi) * P_lm(l, m, theta);
  else if (m < 0)
    return sqrt(2.0) * K_lm(l, m) * sin(abs(m) * phi) * P_lm(l, m, theta);
  else if (m == 0)
    return K_lm(l, m) * P_lm(l, m, theta);

  else
    return -2;
}
