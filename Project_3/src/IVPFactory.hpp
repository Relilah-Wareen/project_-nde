#ifndef IVP_FACTORY_HPP
#define IVP_FACTORY_HPP

#include "ODESolver.hpp"
#include <map>
#include <string>
#include <functional>
#include <memory>

class IVPFactory {
public:
    using CreateIVPCallback = std::unique_ptr<ODESolver> (*)(int, double, std::function<Vector(const Vector&, double)>);

private:
    using CallbackMap = std::map<std::string, CreateIVPCallback>;
    CallbackMap callbacks_;

    IVPFactory() = default;
    IVPFactory(const IVPFactory&) = delete;
    IVPFactory& operator=(const IVPFactory&) = delete;

public:
    static IVPFactory& getInstance() {
        static IVPFactory instance;
        return instance;
    }

    void registerSolver(const std::string& name, CreateIVPCallback createFn) {
        callbacks_[name] = createFn;
    }

    std::unique_ptr<ODESolver> createSolver(
        const std::string& name,
        int p,
        double T,
        std::function<Vector(const Vector&, double)> f) 
    {
        if (callbacks_.find(name) == callbacks_.end()) {
            throw std::runtime_error("Unknown solver type: " + name);
        }
        return callbacks_[name](p, T, f);
    }
};

void RegisterAllSolvers();

#endif // IVP_FACTORY_HPP