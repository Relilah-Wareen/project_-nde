#include "IVPFactory.hpp"
#include "LinearMultistepMethod.hpp"
#include "RKMethod.hpp"
#include "EmbeddedRKMethod.hpp"


namespace {
    std::unique_ptr<ODESolver> createAdamsBashforth(int p, double T, std::function<Vector(const Vector&, double)> f) {
        return std::make_unique<AdamsBashforth>(p, T, f);
    }

    std::unique_ptr<ODESolver> createAdamsMoulton(int p, double T, std::function<Vector(const Vector&, double)> f) {
        return std::make_unique<AdamsMoulton>(p, T, f);
    }

    std::unique_ptr<ODESolver> createBDF(int p, double T, std::function<Vector(const Vector&, double)> f) {
        return std::make_unique<BDFMethod>(p, T, f);
    }

    std::unique_ptr<ODESolver> createFehlberg45(int p, double T, std::function<Vector(const Vector&, double)> f) {
        return std::make_unique<Fehlberg45>(T, f);
    }

    std::unique_ptr<ODESolver> createDormandPrince54(int p, double T, std::function<Vector(const Vector&, double)> f) {
        return std::make_unique<DormandPrince54>(T, f);
    }

    std::unique_ptr<ODESolver> createClassicalRK4(int p, double T, std::function<Vector(const Vector&, double)> f) {
        if (p != 4) throw std::invalid_argument("ClassicalRK4 requires p=4");
        return std::make_unique<ClassicalRK4>(p, T, f);
    }

    std::unique_ptr<ODESolver> createGaussLegendre(int p, double T, std::function<Vector(const Vector&, double)> f) {
        int s = p / 2; 
        if (s < 2 || s > 5) throw std::invalid_argument("GaussLegendre requires p=4,6,8,10");
        return std::make_unique<GaussLegendre>(s, T, f);
    }

    std::unique_ptr<ODESolver> createESDIRK64(int p, double T, std::function<Vector(const Vector&, double)> f) {
        if (p != 4) throw std::invalid_argument("ESDIRK64 requires p=4");
        return std::make_unique<ESDIRK64>(p, T, f);
    }

    std::unique_ptr<ODESolver> createEuler(int p, double T, std::function<Vector(const Vector&, double)> f) {
        if (p != 1) throw std::invalid_argument("Euler requires p=1");
        return std::make_unique<AdamsBashforth>(1, T, f);
    }
}

void RegisterAllSolvers() {
    auto& factory = IVPFactory::getInstance();


    factory.registerSolver("AdamsBashforth", &createAdamsBashforth);
    factory.registerSolver("AdamsMoulton", &createAdamsMoulton);
    factory.registerSolver("BDF", &createBDF);
    factory.registerSolver("Fehlberg45", &createFehlberg45);
    factory.registerSolver("DormandPrince54", &createDormandPrince54);
    factory.registerSolver("ClassicalRK4", &createClassicalRK4);
    factory.registerSolver("GaussLegendre", &createGaussLegendre);
    factory.registerSolver("ESDIRK64", &createESDIRK64);
    factory.registerSolver("Euler", &createEuler);
}
