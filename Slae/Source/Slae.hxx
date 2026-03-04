#include <vector>
#include <stdexcept>
#include <ostream>
class Matrix{
    std::vector<double> mtx_;
public:
    std::size_t nx_;
    std::size_t ny_;
    Matrix(std::size_t nx, std::size_t ny, const std::vector<double>& mtx) : nx_(nx), ny_(ny){
        if (mtx.empty()){
            throw std::invalid_argument("Matrix is empty");
        }
        mtx_ = mtx;
    }
    double operator()(std::size_t i, std::size_t j) const {return mtx_[i * nx_ + j];}

    double& operator()(std::size_t i, std::size_t j) {return mtx_[i * nx_ + j];}

    const std::vector<double>& data() const {return mtx_;}

    Matrix operator+(const Matrix& o) const {
        if(nx_ != o.nx_ or ny_ != o.ny_){throw std::invalid_argument("Sizes of matrix are different");}
        std::vector<double> res = mtx_;
        for(std::size_t i = 0; i < mtx_.size(); ++i){
            res[i] = mtx_[i] + o.mtx_[i];
        }
        return {nx_, ny_, res};
    }
    Matrix operator*(const Matrix& o) const {
        if(nx_ != o.ny_) throw std::invalid_argument("Impossible to multiply");
        std::size_t s = o.nx_ * ny_;
        std::vector<double> res(s, 0);
        for(std::size_t i = 0; i < ny_; ++i){
            for(std::size_t j = 0; j < o.nx_; ++j){
                double current = 0;
                for(std::size_t r = 0; r < nx_; ++r){
                    current += mtx_[i * nx_ + r] * o.mtx_[r * o.nx_ + j];
                }
                res[i * o.nx_ + j] = current;
            }
        }
        return {o.nx_, ny_, res};
    }
    Matrix operator*(double a) const {
        std::vector<double> res = mtx_;
        for(auto& item : res) {item *= a;}
        return {nx_, ny_, res};
    }
};
inline std::ostream& operator<<(std::ostream& os, const Matrix& mtx){
    for(std::size_t i = 0; i < mtx.ny_; ++i){
        for(std::size_t j = 0; j < mtx.nx_; ++j){
            os << mtx(i, j) << " ";
        }
        os << std::endl;
    }
    return os;
}