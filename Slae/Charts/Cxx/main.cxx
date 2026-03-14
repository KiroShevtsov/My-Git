#include <Slae.hxx>
#include <chrono>
#include <fstream>
#include <random>
#include <iostream>
template <typename T> double IterationTime(T&& function, std::size_t iter){
    auto start = std::chrono::steady_clock::now();
    std::forward<T>(function)();
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> tau = end - start;
    return tau.count() / iter;
}
std::vector<double> Mtx(std::size_t nx, std::size_t ny, double alpha, double min = 0, double max = 10){
    std::size_t size = nx * ny;
    std::vector<double> res(size);
    // res.reserve(size);

    std::random_device dev;
    std::mt19937 r(dev());
    std::uniform_real_distribution<double> values(min, max);
    std::uniform_real_distribution<double> prob_non_zero(0, 1);

    for(std::size_t i = 0; i < size; ++i){
        double l = prob_non_zero(r);
        double v = (l <= alpha) ? values(r) : 0;
        res[i] = v;
    }
    return res;
}
int main(){
    try{
        std::size_t nx = 1000;
        std::size_t ny = 1000;

        std::vector<double> k;
        k.reserve(nx);
        for(std::size_t i = 0; i < nx; ++i){
            k.push_back(i * 0.1);
        }
        Vector ksi(k);

        std::ofstream ofs("data.txt");
        if(!ofs.is_open()){
            throw std::runtime_error("didnt open");
        }

        std::vector<double> density;
        std::size_t iter = 100;
        for(std::size_t i = 1; i < iter; ++i){
            density.push_back(i * 1e-2);
        }

        using key = std::pair<std::size_t, std::size_t>;
        
        for(const auto& alpha : density){
            std::vector<double> current = Mtx(nx, ny, alpha);
            Matrix first(nx, ny, current);

            std::map<key, double> m;
            for(std::size_t i = 0; i < ny; ++i){
                for(std::size_t j = 0; j < nx; ++j){
                    double v = current[i * nx + j];
                    if(v != 0){
                        m[std::make_pair(i, j)] = v;
                    }
                }
            }
            SparseMatrix second(nx, ny, m);
            
            double t_dense = IterationTime([&](){Vector r = first * ksi;}, density.size());
            double t_sparse = IterationTime([&](){Vector r = second * ksi;}, density.size());
            ofs << alpha << " " << t_dense << " " << t_sparse << std::endl;
        }
        ofs.close();
    } catch(std::exception& e){
        std::cout << e.what() << std::endl;
    }
}
