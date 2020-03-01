#!/usr/bin/env python3

from scipy import stats, optimize
import numpy as np
import matplotlib.pyplot as plt
import sys

'''Usage: ./approx_beta.py 5 2 100000 0.1 1  0.1'''

class minimize_tvd():
    def __init__(self, beta_a, beta_b, x_lim0, x_lim1, step):
        self.beta_a = beta_a
        self.beta_b = beta_b
        self.x_lim0 = x_lim0
        self.x_lim1 = x_lim1
        self.step = step
        self.x_s = np.arange(self.x_lim0, self.x_lim1, step)
        self.beta_densities = stats.beta.pdf(self.x_s, self.beta_a, self.beta_b)

    def disc_dist(self, phi, mu0, mu1, sigma):
        denom = 1 + ((1/phi) - 1)*np.exp(1/np.square(sigma)*((-1*self.x_s*(mu1-mu0))+(1/2*(np.square(mu1)-np.square(mu0)))))
        return 1/denom

    def tvd(self, density):#total variance distance
        tvd = 0.5*np.sum(np.abs(self.beta_densities-density))
        return tvd

    def objective_function(self, x):
        phi, mu0, mu1, sigma = x[0], x[1], x[2], x[3]
        return(self.tvd(self.disc_dist(phi, mu0, mu1, sigma)))



    def run_min(self, init_phi, init_mu0, init_mu1, init_sigma):
        self.init_phi = init_phi
        self.init_mu0 = init_mu0
        self.init_mu1 = init_mu1
        self.init_sigma = init_sigma
        initial_guesses = np.array([self.init_phi, self.init_mu0, self.init_mu1, self.init_sigma])
        def cb( xk ):
            print( "Current params", xk )
        optimize_result = optimize.minimize(self.objective_function, initial_guesses, method='Nelder-Mead', callback=cb, options=dict(maxiter=10000, maxfev=10000))
        success = optimize_result['success']
        fit_x = optimize_result['x']
        if success == True:
            self.phi_f, self.mu0_f, self.mu1_f, self.sigma_f = fit_x[0], fit_x[1], fit_x[2], fit_x[3]
        else:
            print("else: {}".format(optimize_result["message"]))
            self.phi_f, self.mu0_f, self.mu1_f, self.sigma_f = self.init_phi, self.init_mu0, self.init_mu1, self.init_sigma

    def get_x(self):
        return self.x_s

    def get_beta_y(self):
        return self.beta_densities

    def get_best_fit_densities(self):
        self.disc_densities = self.disc_dist(self.phi_f, self.mu0_f, self.mu1_f, self.sigma_f)
        return self.disc_densities

    def get_best_fit_params(self):
        return (self.phi_f, self.mu0_f, self.mu1_f, self.sigma_f)

    def get_best_tvd(self):
        return self.tvd(self.disc_dist(self.phi_f, self.mu0_f, self.mu1_f, self.sigma_f))


#beta_a = 2.5
#beta_b = 1
# beta_a = 5
# beta_b = 2
beta_a = float(sys.argv[1]) #5
beta_b = float(sys.argv[2]) #2
num_points = int(sys.argv[3]) #100000
x_lim0 = float(sys.argv[4]) #0.01
x_lim1 = float(sys.argv[5]) #1
step = float(sys.argv[6]) #0.01

colors = ['dodgerblue', 'coral', 'darkslategray', 'purple', 'darkorange', 'olive', 'rosybrown','black', 'khaki','grey', 'saddlebrown', 'tan', 'darkgoldenrod', 'turquoise', 'mediumvioletred', 'cyan', 'navy', 'indigo', 'magenta', 'deeppink']
counter = 0

min_tvd = 100000
approx_class = minimize_tvd(beta_a, beta_b, x_lim0, x_lim1, step)
x_s = approx_class.get_x()
beta_y_s = approx_class.get_beta_y()
fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
ax1.plot(x_s, beta_y_s, color='black')
ax3.plot(x_s, beta_y_s, color='black')
for init_phi in stats.uniform.rvs(0.5, 0.5, num_points):
    for init_mu0 in stats.uniform.rvs(-10,10, num_points):
        for init_mu1 in stats.uniform.rvs(0.5, 10, num_points):
            for init_sigma in stats.uniform.rvs(0.25, 10, num_points):
                #try scipy minimize
                approx_class.run_min(init_phi, init_mu0, init_mu1, init_sigma)
                fit_tvd = approx_class.get_best_tvd()
                if fit_tvd < min_tvd:
                    best_fit_y = approx_class.get_best_fit_densities()
                    best_fit_phi, best_fit_mu0, best_fit_mu1, best_fit_sigma = approx_class.get_best_fit_params()
                    min_tvd = fit_tvd
                    ax1.plot(x_s, best_fit_y, color=colors[counter%len(colors)], linestyle='--')
                    ax2.scatter(counter, fit_tvd)
                    counter += 1


print("Phi: {}\nmu0: {}\nmu1: {}\nsigma: {}\nTVD: {}".format(best_fit_phi, best_fit_mu0, best_fit_mu1, best_fit_sigma, min_tvd))
#plot fit
ax3.plot(x_s, best_fit_y, color='darkorange', linestyle='--')
ax3.annotate("Phi: {}\nmu0: {}\nmu1: {}\nsigma: {}\nTVD: {}".format(best_fit_phi, best_fit_mu0, best_fit_mu1, best_fit_sigma, min_tvd), xy=(0.2,0.2))
fig.savefig('disc_opt.png')
plt.close(fig)
