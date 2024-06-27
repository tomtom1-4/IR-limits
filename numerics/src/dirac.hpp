#ifndef DIRAC_HEADER
#define DIRAC_HEADER

static std::vector<std::vector<double> > metric = {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}};

static int LeviCevita(int alpha, int beta, int gamma, int delta) {
    int indices[4] = {alpha, beta, gamma, delta};
    int counter = 0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++){
            if (indices[j] == i and i != j) {
                counter ++;
                indices[j] = indices[i];
                indices[i] = i;
                //std::cout << counter<< ": " << indices[0] << ", " << indices[1] << ", " << indices[2] << ", " << indices[3] << std::endl;
            }
        }
    }
    if (counter % 2 == 0) {
        return -1;
    }
    else {
        return 1;
    }
}

static std::vector<std::vector<std::complex<double> > > gamma0 = {{0, 0, 1, 0},
                                                            {0, 0, 0, 1},
                                                            {1, 0, 0, 0},
                                                            {0, 1, 0, 0}};
static std::vector<std::vector<std::complex<double> > > gamma1 = {{0, 0, 0, 1},
                                                            {0, 0, 1, 0},
                                                            {0, -1, 0, 0},
                                                            {-1, 0, 0, 0}};
static std::vector<std::vector<std::complex<double> > > gamma2 = {{0, 0, 0, -I},
                                                            {0, 0, I, 0},
                                                            {0, I, 0, 0},
                                                            {-I, 0, 0, 0}};
static std::vector<std::vector<std::complex<double> > > gamma3 = {{0, 0, 1, 0},
                                                            {0, 0, 0, -1},
                                                            {-1, 0, 0, 0},
                                                            {0, 1, 0, 0}};

static std::vector<std::vector<std::vector<std::complex<double> > > > gammas = {gamma0, gamma1, gamma2, gamma3};

static std::vector<std::complex<double>> spinor(double *p, int pol) {
    double s = std::abs(p[0])/p[0];
    std::vector<std::complex<double>> u;
    double cos = p[1]/sqrt(p[1] * p[1] + p[2] * p[2]);
    double sin_half = -std::sqrt((1 - cos)/2);
    double cos_half = std::sqrt(1 - sin_half * sin_half);
    std::complex<double> phase = cos_half + I * sin_half;
    phase = 1;
    if (std::abs(p[0] + p[3]) > 1e-12) {
        if (pol == 1) {
            u.push_back(0);
            u.push_back(0);
            u.push_back(sqrt(s * (p[0] + p[3])) * phase);
            u.push_back(s * (p[1] + I * p[2])/sqrt(s * (p[0] + p[3])) * phase);
        }
        else if (pol == -1) {
            u.push_back(s * (-p[1] + I * p[2])/sqrt(s * (p[0] + p[3])) * std::conj(phase));
            u.push_back(sqrt(s * (p[0] + p[3])) * std::conj(phase));
            u.push_back(0);
            u.push_back(0);
        }
        return u;
    }
    else {
        // Failsafe solution for p[0] + p[3] -> 0
        if (pol == 1) {
            u.push_back(0);
            u.push_back(0);
            u.push_back(0);
            u.push_back(sqrt(2 * s * p[0]));
        }
        else if (pol == -1) {
            u.push_back(-sqrt(2 * s * p[0]));
            u.push_back(0);
            u.push_back(0);
            u.push_back(0);
        }
        return u;
    }
}

static std::vector<std::complex<double>> polarization (double *p, int pol) {
    std::vector<std::complex<double>> ep;

    double pT = std::sqrt(p[1] * p[1] + p[2] * p[2]);
    if (std::abs(p[1]) > 1e-12 or std::abs(p[2]) > 1e-12) {
        ep.push_back(0);
        ep.push_back((-pol * p[1] * p[3] - I * p[2] * std::abs(p[0]))/std::sqrt(2)/std::abs(p[0])/pT);
        ep.push_back((-pol * p[2] * p[3] + I * p[1] * std::abs(p[0]))/std::sqrt(2)/std::abs(p[0])/pT);
        ep.push_back(pol * pT/std::abs(p[0])/std::sqrt(2));
        return ep;
    }
    else {
        // Failsafe solution for p[1], p[2] -> 0
        std::cout << "Failsafe solution" << std::endl;

        ep.push_back(0);
        ep.push_back(-pol/sqrt(2));
        ep.push_back(-I/sqrt(2));
        ep.push_back(0);
        return ep;
    }
}


static std::vector<std::complex<double>> complex_conjugate(std::vector<std::complex<double>> u) {
    std::vector<std::complex<double>> out;
    for (int i = 0; i < 4; i++) {
        out.push_back(std::conj(u[i]));
    }
    return out;
}

static std::complex<double> operator*(std::vector<std::complex<double>> u, std::vector<std::complex<double>> v) {
    std::complex<double> out = 0;
    for (int dummy = 0; dummy < 4; dummy++){
        out = out + std::conj(u[dummy]) * v[dummy];
    }
    return out;
}

static std::vector<std::complex<double>> operator+(std::vector<std::complex<double>> u, std::vector<std::complex<double>> v) {
    std::vector<std::complex<double>> out = {0., 0., 0., 0.};
    for (int i = 0; i < 4; i++){
        out[i] = u[i] + v[i];
    }
    return out;
}

static std::vector<std::complex<double> > operator*(std::vector<std::vector<std::complex<double>>> gamma_mu, std::vector<std::complex<double>> u) {
    std::vector<std::complex<double> > out = {0, 0, 0, 0};
    for (int i = 0; i < 4; i++){
        for (int dummy = 0; dummy < 4; dummy++){
            out[i] = out[i] + gamma_mu[i][dummy] * u[dummy];
        }
    }
    return out;
}

static std::vector<std::vector<std::complex<double> > > operator*(std::vector<std::vector<std::complex<double> > > gamma_mu,
                                                         std::vector<std::vector<std::complex<double> > > gamma_nu) {
    std::vector<std::vector<std::complex<double> > > out = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++) {
            for (int dummy = 0; dummy < 4; dummy++) {
                out[i][j] = out[i][j] + gamma_mu[i][dummy] * gamma_nu[dummy][j];
            }
        }
    }
    return out;
}

static std::vector<std::vector<std::complex<double> > > operator+(std::vector<std::vector<std::complex<double> > > gamma_mu,
                                                         std::vector<std::vector<std::complex<double> > > gamma_nu) {
    std::vector<std::vector<std::complex<double> > > out = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++) {
            out[i][j] = gamma_mu[i][j] + gamma_nu[i][j];
        }
    }
    return out;
}

static std::vector<std::vector<std::complex<double> > > operator-(std::vector<std::vector<std::complex<double> > > gamma_mu,
                                                         std::vector<std::vector<std::complex<double> > > gamma_nu) {
    std::vector<std::vector<std::complex<double> > > out = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++) {
            out[i][j] = gamma_mu[i][j] - gamma_nu[i][j];
        }
    }
    return out;
}

template <class A>
static std::vector<std::vector<std::complex<double> > > operator*(std::vector<std::vector<std::complex<double> > > gamma_mu, A a) {
    std::vector<std::vector<std::complex<double> > > out = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            out[i][j] = gamma_mu[i][j] * a;
        }
    }
    return out;
}

static std::vector<std::complex<double> >  operator*(std::vector<std::complex<double> >  p, std::complex<double> a) {
    std::vector<std::complex<double> > out = {0,0,0,0};
    for (int i = 0; i < 4; i++) {
        out[i] = p[i] * a;
    }
    return out;
}

static std::vector<std::complex<double> >  operator*(std::complex<double> a, std::vector<std::complex<double> >  p) {
    std::vector<std::complex<double> > out = {0,0,0,0};
    for (int i = 0; i < 4; i++) {
        out[i] = p[i] * a;
    }
    return out;
}

template <class A>
static std::vector<std::vector<std::complex<double> > > p_slash(A *p) {
    std::vector<std::vector<std::complex<double> > > out = {{0, 0, 0, 0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
    for (int mu = 0; mu < 4; mu++) {
        out = out + gammas[mu] * p[mu] * metric[mu][mu];
    }
    return out;
}

#endif