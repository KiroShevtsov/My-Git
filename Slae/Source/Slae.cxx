#include "Slae.hxx"
/*a - under_diagonal, b - diagonal, c - upper_diagonal for matrix T, where Tx = d*/
Vector SweepMethod(const std::vector<double>& a, const std::vector<double>& b,
    const std::vector<double>& c, const Vector& d){
        std::size_t n = d.dim();
        
        if(b.size() != n or a.size() != (n - 1) or c.size() != (n - 1)){
            throw std::invalid_argument("Vectors don't match tridiagonal matrix elements");
        }
        std::vector<double> alpha(n - 1);
        std::vector<double> beta(n);

        alpha[0] = -(c[0] / b[0]);
        beta[0] = d[0] / b[0];
        for(std::size_t i = 1; i < (n - 1); ++i){
            alpha[i]  = -(c[i]) / (a[i - 1] * alpha[i - 1] + b[i]);
            beta[i] = (d[i] - beta[i - 1] * a[i - 1]) / (a[i - 1] * alpha[i - 1] + b[i]);
        }
        double q = (a[n - 2] * alpha[n - 2] + b[n - 1]);
        beta[n - 1] = (d[n - 1] - a[n - 2] * beta[n - 2]) / q;

        std::vector<double> x(n);
        x[n - 1] = beta[n - 1];
        for(int i = (n - 2); i >= 0; i--){
            x[i] = alpha[i] * x[i + 1] + beta[i];
        }
        return Vector(x);
}